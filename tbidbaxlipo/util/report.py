"""
Code to generate reports documenting model behavior and output.
"""

from tempfile import NamedTemporaryFile
import txt2pdf
from pyPdf import PdfFileWriter, PdfFileReader
import os
from pylab import savefig
from matplotlib.backends.backend_pdf import PdfPages

class Report():
    """A class to facilitate the creation of model analysis documents
    incorporating both figures (from simulations, etc.) and text.
    """
    tempFileList = []

    def addText(self, text):
        """Add text to the report."""
        # Write the text to a temporary text file
        text_file = NamedTemporaryFile(suffix='.txt', mode='w')
        text_file.write(text)
        text_file.flush()
        print text_file.name
        self.tempFileList.append(text_file)

    def addTextFile(self, filename):
        """Include the contents of a text file in the report."""
        text_file = open(filename, 'r')
        text = text_file.read()
        addText(text)
        text_file.close()

    def addCurrentFigure(self):
        """Add the current figure to the report.
        
        Uses the matplotlib function savefig to export the figure to a
        temporary PDF file which is incorporated into the report.
        """

        fig_file = NamedTemporaryFile(suffix='.pdf')
        savefig(fig_file.name)
        #pdf_file = PdfPages(fig_file.name)
        #pdf_file.savefig()
        #pdf_file.close()
        self.tempFileList.append(fig_file)

    def writeReport(self, outputfilename='report'):
        """Output the complete report to the given PDF file."""

        output = PdfFileWriter()

        for tempFile in self.tempFileList:
            tempFileName = tempFile.name
            fileType =  tempFileName[(len(tempFileName)-3):len(tempFileName)]
            input = None

            if fileType == 'pdf':
                input = PdfFileReader(file(tempFileName, "rb"))
            elif fileType == 'txt':
                pdftext = txt2pdf.pyText2Pdf()
                pdftext._ifile = tempFileName
                pdftext.Convert()
                input = PdfFileReader(file(pdftext._ofile, "rb"))
                os.unlink(pdftext._ofile)

            if input:
                TotPgNum = input.getNumPages()

                for i in range(TotPgNum):
                    output.addPage(input.getPage(i))

        outputStream = file(outputfilename + '.pdf', "wb")
        output.write(outputStream)
        outputStream.close()

