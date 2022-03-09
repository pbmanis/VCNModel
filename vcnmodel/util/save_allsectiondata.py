#!/usr/bin/env python
# coding: utf-8
# originally from an ipynb.
# In[1]:


# get_ipython().run_line_magic('matplotlib', 'qt')
import matplotlib
from dataclasses import dataclass, field
import re
from typing import Union, List, Tuple
import toml

import datetime as datetime
from dateutil.tz import tzlocal
import numpy as np

matplotlib.use("Qt5Agg")
from pathlib import Path

import pandas as pd
import pynwb as NWB
from matplotlib import pyplot as mpl

config = toml.load(open("wheres_my_data.toml", "r"))

re_secn = re.compile("create sections\[(?P<nsections>\d+)\]")
re_swcs = re.compile("// SWC ID = (?P<swcid>\d+)")


def defemptylist():
    return []


def defcliplist():
    return [180., 300.0]


@dataclass
class ModelData:
    clip_window: List = field(
        default_factory=defcliplist
    )  # clipping window for data time.
    time: List = field(default_factory=defemptylist)  # 1d np array for time
    data: List = field(default_factory=defemptylist)  # 2d np array section x time
    sections: List = field(default_factory=defemptylist)  # 1d np array section (string)
    cmd_wave: List = field(
        default_factory=defemptylist
    )  # for IV (2d np array section x time like data)
    an_spikes: List = field(default_factory=defemptylist)
    params: str = ""


"""
Plot all locations (for reference) for one stimulus (AN data)
Each trial contains a dict:
dict_keys(['stimInfo', 'spikeTimes', 'inputSpikeTimes', 'time',
'somaVoltage', 'dendriteVoltage', 'allDendriteVoltages', 
'stimWaveform', 'stimTimebase'])

"""


def plot_all_locs_AN(MD):
    time = MD.time
    data_array = MD.data
    f, ax = mpl.subplots(1, 1)
    spines_to_remove = ["top", "right"]
    ax = [ax]
    for a in ax:
        for spine in spines_to_remove:
            a.spines[spine].set_visible(False)
    # for trial in d['trials']:
    # print(len(data_array))
    for i in range(len(data_array)):
        ax[0].plot(time, data_array[i, :], linewidth=0.25)

    mpl.show()


"""
Plot all locations for one stimulus (IV data)
"""


def plot_all_locs_IV():
    # f, ax = mpl.subplots(1,1)
    spines_to_remove = ["top", "right"]
    ax = [ax]

    for a in ax:
        for spine in spines_to_remove:
            a.spines[spine].set_visible(False)
    for i, R in enumerate(d["Results"]):
        if i < len(d["Results"]) - 1:
            continue
        #     print(R)
        inj = list(R.keys())[0]
        r = R[inj]
        for v in r["vec"]:
            ax[0].plot(r["monitor"]["time"], r["vec"][v])
    mpl.show()



# """
# Save IV data to csv file with time on left column and sections ordered to the right
# we do this by building the csv as a pandas data frame and filling it
# """
#
#
# def save_to_csv():
#     for i, R in enumerate(d["Results"]):
#         if i < len(d["Results"]) - 1:
#             continue
#         #     print(R)
#         inj = list(R.keys())[0]
#         r = R[inj]
#         secnos = len(r["vec"])
#         secs = ["time"]
#         for s in range(secnos):
#             secs.append(f"sections[{s:d}]")
#         df = pd.DataFrame(columns=secs)
#         print(df.columns)
#         df["time"] = r["monitor"]["time"]
#         for v in secs[1:]:
#             df[v] = r["vec"][v]
# df.to_csv(Path('..', fn.name).with_suffix('.csv'))
# df.to_json(Path('..', fn.name).with_suffix('.json'))
# df.to_hdf(Path('..', fn.name).with_suffix('.hdf'), key='Voltage')
# df.to_parquet(Path('..', fn.name).with_suffix('.parquet'))


"""
Save AN data to NWB formatted file

Another (and as I think about it, probably preferred approach) would be 
to do this by adding icephys_electrodes and populating it with a list of
IntracellularElectrodes (icephys.IntracellularElectrode) that use 
the “location” str field to store the map back to the SWC for each HOC section. 
Although this will take a little more space, it makes more sense in terms
of the actual experiment/model, and could potentially be either sparse 
or fully populated. The structure appears to be able to hold multiple 
electrodes in a list (at least the documentation leads me to believe that),
so it is still just one read. I say it is logical because the simulation
involves inserting electrodes into each section, as if this were a real cell
and we could actually get that many electrodes into it without damage.

"""


def save_AN_to_NWB(
    cell,
    filename: Union[Path, str],
    modeldata: object,
    params: dict,
    swcmap=dict,
    outpath: Union[Path, None] = None,
    sessionid: int = 0,
):

    outfilename = "test"
    datafilename = filename
    sessionno = 0
    subject = NWB.file.Subject(
        age=str(0),
        description="vcnmodel",
        genotype="None",
        sex="None",
        species="Computer",
        subject_id=f"VCN_c{cell:02d}",
        weight="None",
    )

    #     print(info)
    nwbfile = NWB.NWBFile(
        "VCN Data set",
        str(datafilename),
        datetime.datetime.now(tzlocal()),
        experimenter="Manis, Paul B.",
        lab="Manis Lab",
        institution="UNC Chapel Hill",
        experiment_description="Model Output",
        session_id=f"{sessionid:d}",
        notes=modeldata.params,
        protocol=f"Cell Type: Bushy XM13_nacncoop AN Model",
        subject=subject,
    )
    device = NWB.device.Device(
        "NEURON",
        manufacturer="Moore, Hines and Carnevale",
        description="NEURON7.7 and Python3.7",
    )
    nwbfile.add_device(device)

    electrode_type = NWB.icephys.IntracellularElectrode(
        name="site electrode",
        device=device,
        description="Neuron section voltage. See icephys_electrode structure for mapping NEURON sections back to SWC reconstruction",
        resistance="0.0",  # no series resistance
    )
    nwbfile.add_icephys_electrode(electrode_type)

    electrodes = [None] * len(modeldata.sections)
    for i, section in enumerate(modeldata.sections):
        # print("Section: ", section)
        electrodes[i] = NWB.icephys.IntracellularElectrode(
            name=section,  # neuron section name
            description="Intracellular electrode",
            device=device,
            location=str(
                swcmap[section]
            ),  # all of the swc elements in the neuron section
        )
    nwbfile.add_icephys_electrode(electrodes)  # add all of the electrode information
    # print(modeldata.params)
    istim = NWB.icephys.CurrentClampStimulusSeries(
        name="Ics",
        data=np.array(modeldata.cmd_wave),
        unit="amperes",
        starting_time=0.0,  # info["__timestamp__"],
        rate=params.dtIC, # ["dt"],
        electrode=electrode_type,
        gain=1.0,
        sweep_number=np.uint64(1),  # othwerwise throws warning in NWB
    )

    vdata = NWB.icephys.CurrentClampSeries(
        name="Vcs",
        data=modeldata.data,
        unit="volts",
        electrode=electrode_type,
        gain=1.0,
        bias_current=0.0,
        bridge_balance=0.0,
        capacitance_compensation=0.0,
        stimulus_description="NA",
        resolution=0.0,
        conversion=1.0,
        timestamps=None,
        starting_time=modeldata.clip_window[0],  # info["__timestamp__"],
        rate=params.dtIC, # ["dt"],
        comments="no comments",
        description="no description",
        control=None,
        control_description=None,
        sweep_number=np.uint64(1),  # othwerwise throws warning in NWB
    )
    ANSpikeOrigin = []
    for i in range(len(modeldata.an_spikes)):
        print(i, modeldata.an_spikes[i])
        ANSpikeOrigin.append(i * np.ones(len(modeldata.an_spikes[i])))
    an_spike_origin = [y for x in ANSpikeOrigin for y in x]
    an_time_stamps = [y for x in modeldata.an_spikes for y in x]
    print('len ts: ', len(an_time_stamps))
    print('anspkorigin: ', an_spike_origin)
    aninputs = NWB.base.TimeSeries(
        name="ANSpikeTimes",
        data=an_spike_origin,
        unit="ms",
        comments="Input AN spike trains",
        description="AN Spike Trains; data is stored in TimeStamps as a dictionary",
        timestamps=an_time_stamps,
    )

    nwbfile.add_stimulus(aninputs)
    nwbfile.add_acquisition(istim)
    nwbfile.add_acquisition(vdata)

    if outpath is not None:
        fout = str(Path(outpath, outfilename)) + ".nwb"
        with NWB.NWBHDF5IO(fout, "w") as io:
            io.write(nwbfile)
            print(f"Ok: {fout:s}")


def read_nwb(infilename):
    """
    Read the nwb file back in and display the data
    Note that the [()] construct is essential for accessing the data

    """

    infile = Path(infilename).with_suffix(".nwb")
    if infile.is_file():
        print("found file: ", infile)
        try:
            io = NWB.NWBHDF5IO(str(infile), "r")
            nwbdata = io.read()
            notes = nwbdata.notes
            protocol = nwbdata.protocol
            vcs = nwbdata.acquisition["Vcs"]
            data = vcs.data[()]
            ANSpikes = nwbdata.stimulus["ANSpikeTimes"]
            electrodes = nwbdata.icephys_electrodes
            anspikeorigin = ANSpikes.data[()]
            anspiketimes = ANSpikes.timestamps[()]
        except:
            print("Error reading data file")
            exit()
        finally:
            io.close()

    for elecname in electrodes:
        elec = electrodes[elecname]

    timebase = np.arange(0, vcs.rate * data.shape[1], vcs.rate)+vcs.starting_time
    xlims = [np.min(timebase), np.max(timebase)]

    f, axl = mpl.subplots(2, 1, sharex=True)
    for i in range(data.shape[0]):  # plot traces for every section
        axl[0].plot(timebase, data[i, :], linewidth=0.35)

    ost = np.argsort(anspiketimes)  # ordered sort times
    anspikeorigin = anspikeorigin[ost]
    anspiketimes = anspiketimes[ost]
    ninputs = np.max(anspikeorigin)
    coincidence_window = 0.05  # adjust this to change the coincidence window
    intervals = np.where(np.diff(anspiketimes) <= coincidence_window)[0]

    # plot AN spikes that occur together within a close coincidence window (like Joris and Smith)
    axl[1].scatter(
        anspiketimes[intervals],
        anspikeorigin[intervals],
        c="r",
        s=25,
        marker="|",
        linewidth=2,
        facecolor="r",
        edgecolor="r",
    )
    axl[1].scatter(
        anspiketimes[intervals + 1],
        anspikeorigin[intervals + 1],
        c="r",
        s=25,
        marker="|",
        linewidth=2,
        facecolor="r",
        edgecolor="r",
    )

    noncoinc = [i 
        for i, x in enumerate(anspiketimes)
        if x not in intervals
    ]

    axl[1].scatter(
        anspiketimes[noncoinc]+1,
        anspikeorigin[noncoinc],
        c="k",
        s=25,
        marker="|",
        linewidth=1,
        facecolor="k",
        edgecolor="None",
    )

    # xlim = axl[0].get_xlim()
    axl[0].set_xlim(xlims)
    axl[1].set_xlim(xlims)
    axl[0].set_ylabel("V(mv)", fontsize=9)
    axl[1].set_ylabel("Input \#", fontsize=9)
    axl[1].set_xlabel("T (ms)")
    axl[0].set_title("All section voltages")
    axl[1].set_title(f"Input Spikes and Coincidences (dt={coincidence_window:.2f} ms)")

    for ax in axl:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    mpl.show()


def write_file():

    cell = 11
    hocfile = f"{config['cellDataDirectory']:s}/VCN_c{cell:02d}/Morphology/VCN_c{cell:02d}_Full.hoc"
    basepath = Path(config['baseDataDirectory'])
    opath = Path(basepath, "ForGeorge")
    with open(hocfile, "r") as fh:
        hocdata = fh.read()
    gr = re.search(re_secn, hocdata)

    nsections = int(gr.group("nsections"))
    swcmap = {}
    for n in range(nsections):
        sec_str = "sections\[{0}\]".format(n)
        re_secdata = re.compile(sec_str + " {(?P<hocs>.*?)}", re.DOTALL)
        line = re.findall(re_secdata, hocdata)
        if line is None:
            swcmap[secstr] = None
            continue
        swcs = []

        for li in line[0].split("\n"):
            swcgr = re.search(re_swcs, li)
            if swcgr is not None:
                swcs.append(swcgr.group("swcid"))
        swcl = ",".join(str(s) for s in swcs)
        swcmap[sec_str.replace("\\", "")] = swcl

    fn = Path(
        f"runANPSTH-all-2020-09-09.11-24-50/AN_Result_VCN_c11_inp=self_XM13A_nacncoop_II_soma=1.203_dend=1.510_normal_delays_multisite_001_tonepip_020dB_16000.0_HS.p"
    )
    basepath = Path(config['cellDataDirectory'],'VCN_c11/Simulations/AN/')
    fnb = Path(basepath, fn)
    d = pd.read_pickle(fnb)

    trial = 0  # for now, just handle one trial
    dx = d['Results'][trial]
    MD = ModelData()
    MD.params = toml.dumps(d["Params"].__dict__)  # convert to a string
    twin = np.nonzero(
        (dx["time"] > MD.clip_window[0]) & (dx["time"] <= MD.clip_window[1])
    )
    MD.time = dx["time"][twin]

    data_array = None

    dxl = dx["allDendriteVoltages"].tolist()
    nsec = len(dxl)# print(dxl)
    for sec in dxl:  # repetition

        y = np.array(dxl[sec])[twin]# [twin]
        if len(MD.data) == 0:
            MD.data = np.zeros((nsec, len(y)))
            MD.sections = [None] * nsec
        secno = int(sec[9:-1])  # get the section number
        MD.data[secno, :] = y  # store data in section order
        MD.sections[secno] = str(sec)
    ANSpikeTimes = []  # store as dict so that input number is an explicit key
    for i in range(len(dx["inputSpikeTimes"])):
        ANSpikeTimes.append(dx["inputSpikeTimes"][i])

    for k in range(len(ANSpikeTimes)):
        indx = np.nonzero(
            (ANSpikeTimes[k] > MD.clip_window[0]) & (ANSpikeTimes[k] <= MD.clip_window[1])
        )
        ANSpikeTimes[k] = ANSpikeTimes[k][indx] # + MD.clip_window[0]
        print(ANSpikeTimes[k])
    MD.an_spikes = ANSpikeTimes
    plot_all_locs_AN(MD)
    save_AN_to_NWB(
        cell,
        filename=fn,
        modeldata=MD,
        params=d["Params"],
        swcmap=swcmap,
        outpath=opath,
    )


def read_file():
    fin = Path(config['baseDataDirectory'], "ForGeorge/test.nwb")
    read_nwb(fin)
    return


if __name__ == "__main__":
    # write_file()
    read_file()
