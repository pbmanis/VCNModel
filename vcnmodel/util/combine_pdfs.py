"""
Merge PDF files

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2020 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

from pathlib import Path

from PyPDF2 import PdfFileMerger, PdfFileReader


def merge_pdfs(datatype, outname=None):
    """
    Merge the PDFs in tempdir with the outputfile (self.outputFilename)
    the PDFs are deleted once the merge is complete.
    """
    fnsg = Path(".").glob(f"{datatype:s}*.pdf")

    pdf_merger = PdfFileMerger()
    fns = sorted(list(fnsg))  # list filenames in the tempdir and sort by name
    rdrs = []
    print(fns)
    for fn in fns:
        with open(str(fn), "rb") as fh:
            rdrs.append(fh)
            pdf_merger.append(PdfFileReader(fh))  # now append the rest
    with open(outname, "wb") as fo:
        pdf_merger.write(fo)
    for r in rdrs:  # close the files... though the merger might also do this
        r.close()


def main():
    for datatype in ["Original", "Inflated"]:
        merge_pdfs(datatype, outname=Path(f"{datatype:s}_PSTH.pdf"))


if __name__ == "__main__":
    main()
