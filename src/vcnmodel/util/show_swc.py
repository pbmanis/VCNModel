"""
Render reconstructions and place them into a 3-D space

"""
import argparse
import sys
from pathlib import Path

import neuronvis.swc_to_hoc
import numpy as np
import pandas as pd
from matplotlib import pyplot as mpl
from mayavi import mlab
from neuronvis import hocRender as HR

import toml
from vcnmodel.util.parse_hoc import ParseHoc

config = toml.load(open("wheres_my_data.toml", "r"))

# getting from new soma reconstructions:
# modes:
#   'cs' : cells + swc
#   'c' : just cells
#   's'  : just swcs
#   'n'  : just pointinfo from new coordinates (axon, dendrite(s) connection points, estimated cell center)
#   'x'  : use translated instead of 'swcs' in hoc (appends _translated to stem filename)
#          requires that you have run with 'ns -t' first to generate the translated files
#   'u'  : convert from xlsx to swc then translate
# combinations ('cs', for example is cells + swcs)


def main():
    parser = argparse.ArgumentParser(description="SWC/HOC renderer")
    parser.add_argument(
        type=str,
        default="n",
        dest="mode",
        # choices=['s', 'c', 'n']
        help="mode",
    )
    parser.add_argument(
        "-n",
        dest="cellnames",
        type=str,
        default="None",
        nargs="+",
        help="Specificy a single cell name when plotting",
    )
    parser.add_argument(
        "-t",
        "--translate",
        action="store_true",
        dest="translate",
        help="compute translation between cells",
    )
    parser.add_argument(
        "-o",
        "--original",
        action="store_true",
        dest="original",
        help="display original hoc on top of new one for a cell",
    )

    parser.add_argument(
        "-f",
        "--file",
        type=str,
        default="None",
        help="just display the specified file in the SWC dataset",
    )

    parser.add_argument(
        "-F",
        "--FILENAME",
        type=str,
        default="None",
        dest="filename",
        help="just display the specified file on the path",
    )
    parser.add_argument(
        "--topology",
        action="store_true",
        dest="topology",
        help="Just print topology via neuron",
    )

    args = parser.parse_args()
    mode = (
        args.mode
    )  # mode contains letters 's' (swc), 'n' (new coords) or 'c' (cell reconstruction)
    showOne = False
    allfiles = []
    renderer = "mayavi"
    pmode = "cylinders"  # cylinders or graph or both
    fax = None

    basedir = config["cellDataDirectory"]
    newcoords = Path(basedir, "VCN_Coordinates_16Cells.xlsx")

    newcoorddir = Path(config["baseDataDirectory"], "VCN_newcoord")

    Path.mkdir(newcoorddir, exist_ok=True)

    coorfile = pd.read_excel(newcoords, sheet_name="Sheet1")

    PH = ParseHoc()
    if args.topology and args.filename is not None:
        import neuron
        from neuron import h

        neuron.h.hoc_stdout(
            "/dev/null"
        )  # prevent junk from printing while reading the file
        success = neuron.h.load_file(str(args.filename))
        neuron.h.hoc_stdout()

        h.topology()
        exit()

    if args.file is not "None":
        f = Path(basedir, args.file)

        renderer = "mayavi"
        fighandle = mlab.figure(bgcolor=(0, 0, 0), fgcolor=(1, 1, 1), size=(800, 800))
        HR.Render(
            display_mode="cylinders",
            display_renderer=renderer,
            hoc_file=f,
            fax=fax,
            fighandle=fighandle,
            somaonly=True,
            color=(0.7, 0.7, 0.7),
            label=args.file,
            flags=["notext", "norefline", "norefaxes"],
        )
        fswc = Path(config["baseMorphologyDirectory"], "ASA", "CellBodySWCs")
        fn = Path(args.file).stem[:7] + "_CellBody01.hoc"
        fswc = Path(fswc, fn)
        HR.Render(
            display_mode="cylinders",
            display_renderer=renderer,
            hoc_file=fswc,
            fax=fax,
            fighandle=fighandle,
            somaonly=True,
            color=(0.9, 0.0, 0.0),
            label="SWC Soma",
            flags=["notext", "norefline", "norefaxes"],
        )

        mlab.show()
        exit()

    if args.filename is not "None":
        f = Path(args.file)

        renderer = "mayavi"
        print("args.filename: ")
        fighandle = mlab.figure(bgcolor=(0, 0, 0), fgcolor=(1, 1, 1), size=(800, 800))
        HR.Render(
            display_mode="cylinders",
            display_renderer=renderer,
            hoc_file=args.filename,
            fax=fax,
            fighandle=fighandle,
            somaonly=True,
            color=(0.7, 0.7, 1.0),
            label=Path(args.filename).name,
            flags=["notext", "norefline", "norefaxes"],
        )
        mlab.show()
        exit()

    def align_cells_soma(f0, f1, reorder=False):
        """
        Compute the linear transform that moves f1 soma to f0 soma (center)
        """
        p1 = PH.get_soma(f1)  # swc files
        p0 = PH.get_soma(f0)  # new coordinate files
        if p0 is None:
            return
        if p1 is None:
            return
        # do all the combinations because we don't know which is corre t
        s1_top = np.array(p1[0][0:3])
        s0_top = np.array(p0[0][0:3])
        s0_bottom = np.array(p0[1][0:3])
        s1_bottom = np.array(p1[1][0:3])
        # print('\np1: ', p1, p0)
        # print(f"swctop: {f1n:s}  {str(s1_top):s}")
        # print(f"newtop: {f0n:s} {str(s0_top):s}")
        # print(f"swcbottom: {f1n:s}  {str(s1_bottom):s}")
        # print(f"newbottom: {f0n:s} {str(s1_bottom):s}")
        # print('s1top-s0top: ', s1_top-s0_top)
        trans_xyz = s1_top - s0_bottom
        print(f0, f1)
        print(
            "s1top-s0bottom: ", s1_top - s0_bottom
        )  # works best for the 5 cells that match so far
        # print('s1bottom-s0top: ', s1_bottom-s0_top)
        # print('s1bottom-s0bottom: ', s1_bottom-s0_bottom)
        PH.translate_cell(trans_xyz, fn=f1, reorder=reorder)

    # print(coorfile.columns.values)
    # construct hoc files
    if "u" in mode:
        morphfiles0 = []
        for i, cell in enumerate(coorfile["Cell Number"].values):
            cellfilename = f"VCN_c{cell:02d}.hoc"
            print(cellfilename)
            ft = """objref undefined
    undefined = new SectionList()
    objref soma
    soma = new SectionList()
    objref axon
    axon = new SectionList()
    objref dendrite
    dendrite = new SectionList()
    objref apical_dendrite
    apical_dendrite = new SectionList()

    create sections[1]
    access sections[0]
    soma.append()"""
            cxs = coorfile["Cell Body center"].values[i]
            cx_data = [float(x.strip()) for x in cxs.split(",")]
            axs = coorfile["Axon base"].values[i]
            ax_data = [float(x.strip()) for x in axs.split(",")]
            dxs = coorfile["Dendrite 1 base"].values[i]
            dx_data = [float(x.strip()) for x in dxs.split(",")]

            ft += """\nsections[0] {"""
            ft += f"\n pt3dadd({ax_data[0]:f}, {ax_data[1]:f}, {ax_data[2]:f}, 9.)"
            ft += f"\n pt3dadd({cx_data[0]:f}, {cx_data[1]:f}, {cx_data[2]:f}, 20.)"
            ft += f"\n pt3dadd({dx_data[0]:f}, {dx_data[1]:f}, {dx_data[2]:f}, 2.)"
            ft += "\n}\n"

            # print(ft)
            newf = Path(newcoorddir, cellfilename)
            with open(newf, "w") as fh:
                fh.write(ft)
            # morphfiles0.append(newf)
        # allfiles.extend(morphfiles0)

    colors = [
        (1, 0, 1),
        (0, 1, 0),
        (1, 1, 0),
        (0, 0, 1),
        (0, 1, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0.6, 0.6, 0),
        (0.0, 0.6, 0.6),
        (0.6, 0, 0.6),
    ] * 8
    colormap = {}

    # first get lists of the source data hoc files
    if "n" in mode:
        morphfiles_nc = list(newcoorddir.glob("*.hoc"))
        allfiles.extend(morphfiles_nc)

    if "s" in mode:  # get SWCs to plot
        dataPath = Path(config["cellMrophologyDirectory"], "ASA", "CellBodySWCs")
        morphfiles_swc = list(dataPath.glob("*.hoc"))
        allfiles.extend(morphfiles_swc)

    if "x" in mode:  # get translated SWCs to plot
        dataPath = Path(
            config["cellMrophologyDirectory"], "VCNModel", "ASA", "CellBodySWCs"
        )
        morphfiles_swc = list(dataPath.glob("*_translated.hoc"))
        allfiles.extend(morphfiles_swc)

    if "c" in mode:  # get Cells and swcs
        if args.cellnames[0] == "test.hoc":
            morphfiles_cell = [Path(args.cellname[0])]
        elif len(args.cellnames) > 1:
            morphfiles_cell = []
            dataPath = Path(config["cellDataDirectory"])
            # getting from original data sets (full cells):
            for na in args.cellnames:
                morphfiles_cell.append(
                    Path(dataPath, f"VCN_c{na:s}", "Morphology", f"VCN_c{na:s}.hoc")
                )
        else:
            dataPath = Path(config["cellDataDirectory"])
            # getting from original data sets (full cells):
            cells = dataPath.glob("VCN_c*")
            cells = sorted([c for c in cells if len(c.name) == 7])  # just base cells
            # print(cells)
            morphfiles_cell = []
            for c in cells:
                mfile = c.glob("Morphology/*.hoc")
                mfiles = list(mfile)
                cn = c.name
                for m in mfiles:
                    # print('m.name: ', m.name)
                    if m.stem == cn:
                        morphfiles_cell.append(m)
        allfiles.extend(morphfiles_cell)

    print("allfiles: ", allfiles)

    if (
        (("s" in mode) or ("x" in mode)) and "c" in mode and len(mode) == 2
    ):  # match colors between the same cells, make unmatched cells faint and either "reddish" or "bluish"
        nf = 0
        for f1 in morphfiles_swc:
            f1n = str(f1.stem)
            # print('f1n: ', f1n)
            for f2 in morphfiles_cell:
                f2n = str(f2.stem)
                # print('f2n: ', f2n)
                if f1n.startswith(f2n):
                    colormap[f1n] = colors[nf]
                    colormap[f2n] = colors[nf]
                    nf += 1
            if f1n not in list(colormap.keys()):
                colormap[f1n] = (0.4, 0.25, 0.25)
        for f2 in morphfiles_cell:
            f2n = str(f2.stem)
            if f2n not in list(colormap.keys()):
                colormap[f2n] = (0.25, 0.25, 0.4)

    elif (("s" in mode) or ("x" in mode)) and len(
        mode
    ) == 1:  # plot swc only ust map colors through cells sequentially
        nf = 0
        for f1 in morphfiles_swc:
            f1n = str(f1.stem)
            # print('f1n: ', f1n)
            colormap[f1n] = colors[nf]
            nf += 1
            PH.get_soma(f1)

    elif "c" in mode and len(mode) == 1:  # just map colors through cells sequentially
        nf = 0
        for f1 in morphfiles_cell:
            f1n = str(f1.stem)
            # print('f1n: ', f1n)
            colormap[f1n] = colors[nf]
            nf += 1

    elif "n" in mode and len(mode) == 1:  # new coords
        nf = 0
        for f1 in morphfiles_nc:
            f1n = str(f1.stem)
            # print('f1n: ', f1n)
            colormap[f1n] = colors[nf]
            nf += 1

    elif (
        "n" in mode and (("s" in mode) or ("x" in mode)) and len(mode) in [2, 3]
    ):  # match colors between simple swc files and the new coordinates in morphfiles0
        nf = 0
        for f1 in morphfiles_swc:
            f1n = str(f1.stem)
            # print('f1n: ', f1n)
            for f0 in morphfiles_nc:
                f0n = str(f0.stem)
                # print('f2n: ', f2n)
                if f1n.startswith(f0n):
                    colormap[f1n] = tuple([c * 0.6 for c in colors[nf]])
                    colormap[f0n] = colors[nf]
                    # compute transform : how to move swc to nc position
                    if args.translate:
                        if "n" in mode:
                            reord = True
                        else:
                            reord = False
                        # reord = False
                        align_cells_soma(f0, f1, reorder=reord)

                    nf += 1
            if f1n not in list(colormap.keys()):
                colormap[f1n] = (0.4, 0.25, 0.25)
        for f0 in morphfiles_nc:
            f0n = str(f0.stem)
            if f0n not in list(colormap.keys()):
                colormap[f0n] = (0.25, 0.25, 0.4)
        if args.translate:
            print("Cell translations are complete; re-run to display")
            exit()
    elif (
        "n" in mode and "c" in mode and len(mode) in [2, 3]
    ):  # match colors between simple swc files and the new coordinates in morphfiles0
        nf = 0
        for f1 in morphfiles_cell:
            f1n = str(f1.stem)
            # print('f1n: ', f1n)
            for f0 in morphfiles_nc:
                f0n = str(f0.stem)
                # print('f2n: ', f2n)
                if f1n.startswith(f0n):
                    colormap[f1n] = colors[nf]
                    colormap[f0n] = colors[nf]
                    nf += 1
            if f1n not in list(colormap.keys()):
                colormap[f1n] = (0.7, 0.7, 0.7)
        for f0 in morphfiles_nc:
            f0n = str(f0.stem)
            if f0n not in list(colormap.keys()):
                colormap[f0n] = (0.4, 0.4, 0.4)

    elif (
        (("s" in mode) or ("x" in mode)) and "c" in mode and len(mode) == 2
    ):  # match colors between simple swc files and the new coordinates in morphfiles0
        nf = 0
        for f1 in morphfiles_swc:
            f1n = str(f1.stem)
            # print('f1n: ', f1n)
            for f0 in morphfiles_cell:
                f0n = str(f0.stem)
                # print('f2n: ', f2n)
                if f1n.startswith(f0n):
                    colormap[f1n] = colors[nf]
                    colormap[f0n] = colors[nf]
                    nf += 1
            if f1n not in list(colormap.keys()):
                colormap[f1n] = (0.7, 0.7, 0.7)
        for f0 in morphfiles_nc:
            f0n = str(f0.stem)
            if f0n not in list(colormap.keys()):
                colormap[f0n] = (0.4, 0.4, 0.4)
    else:
        raise ValueError("mode: %s not implemented", mode)

    # print([f for f in swcFiles])

    if renderer.startswith("mpl"):
        fig = mpl.figure()
        ax = fig.add_subplot(111, projection="3d")
        fax = [fig, ax]
        fighandle = fig
        renderer = "mpl"
    elif renderer.startswith("mayavi"):
        renderer = "mayavi"
        fighandle = mlab.figure(bgcolor=(0, 0, 0), fgcolor=(1, 1, 1), size=(800, 800))

    else:
        renderer = "pyqtgraph"
        fighandle = None

    renders = []
    print(allfiles)
    for i, f in enumerate(allfiles):
        fn = str(f.stem)
        print("FN: ", fn)
        # if args.cellnames is not 'None':
        #     for i, f in enumerate(args.cellnames):
        #         print(fn, args.cellnames[i])
        #     if not args.cellnames[i].startswith('test') and not fn.endsswith(args.cellnames[i]):
        #         print('not available')
        #         continue
        if "_untranslated" in fn:
            continue
        gmode = pmode
        if pmode != "both":
            usemode = pmode
            if f.parent.match("*VCN_newcoord"):  # just draw the lines
                usemode = "graph"
                colormap[fn] = (
                    1,
                    1,
                    1,
                )  # and do them in white - we know they are already in position (mostly)
            renders.append(
                HR.Render(
                    display_mode=usemode,
                    display_renderer=renderer,
                    hoc_file=f,
                    fax=fax,
                    fighandle=fighandle,
                    somaonly=True,
                    color=colormap[fn],
                    label=fn,
                    flags=["noaxes", "lines", "notext"],
                )
            )
            if args.original:
                fp = f.parent
                fs = f.stem
                forig = Path(fp, fs + "_original.hoc")
                renders.append(
                    HR.Render(
                        display_mode=usemode,
                        display_renderer=renderer,
                        hoc_file=forig,
                        fax=fax,
                        fighandle=fighandle,
                        somaonly=True,
                        color=colormap[fn],
                        label=fn,
                        flags=["noaxes", "nolines", "notext"],
                    )
                )

        else:
            renders.append(
                HR.Render(
                    display_mode="cylinders",
                    display_renderer=renderer,
                    hoc_file=f,
                    fax=fax,
                    fighandle=fighandle,
                    somaonly=True,
                    color=colormap[fn],
                    label=fn,
                    flags=["noaxes", "nolines", "notext"],
                )
            )
            renders.append(
                HR.Render(
                    display_mode="graph",
                    display_renderer=renderer,
                    hoc_file=f,
                    fax=fax,
                    fighandle=fighandle,
                    somaonly=True,
                    color=(1.0, 1.0, 1.0),
                    label=fn,
                    flags=["noaxes", "nolines", "notext"],
                )
            )
        if i == 0 and renderer == "pyqtgraph":
            fighandle = renders[-1].view
        if i == 0 and showOne:
            break
        print("Plotted: ", fn)

    if renderer == "mpl":
        mpl.show()
    if renderer == "mayavi":
        mlab.show()
    else:
        if sys.flags.interactive == 0:
            import pyqtgraph as pg

            pg.Qt.QtGui.QApplication.exec_()


if __name__ == "__main__":
    main()
