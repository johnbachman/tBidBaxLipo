"""
Code to generate reports documenting model behavior and output.
"""

from tempfile import NamedTemporaryFile
import txt2pdf
from pyPdf import PdfFileWriter, PdfFileReader
import os
import os.path
from pylab import savefig
from matplotlib.backends.backend_pdf import PdfPages
import sys
from pygments import highlight
from pygments.lexers import PythonLexer
from pygments.formatters import LatexFormatter
import re
import subprocess

class Report():
    """A class to facilitate the creation of model analysis documents
    incorporating both figures (from simulations, etc.) and text.
    """

    def __init__(self):
        self.tempFileList = []

    def add_python_code(self, filename):
        """Adds the given Python file to the report after applying
        syntax highlighting with pygments."""

        input_file = open(filename, 'r')
        tex_file = NamedTemporaryFile(suffix='.tex', mode='w')
        tex_file_name = os.path.basename(tex_file.name)
        pdf_filename = re.sub('\.tex$', '.pdf', tex_file_name)

        code = input_file.readlines()
        code = ''.join(code)
        formatter = LatexFormatter(style='colorful', full=True)

        highlight(code, PythonLexer(), formatter, outfile=tex_file)
        tex_file.flush()

        #print 'fname ' + tex_file.name
        p = subprocess.Popen(['pdflatex', tex_file.name],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (p_out, p_err) = p.communicate()
        if p.returncode:
            raise Exception(p_out.rstrip("at line")+"\n"+p_err.rstrip())

        self.tempFileList.append(open(pdf_filename, 'r'))

        os.unlink(re.sub('\.pdf$', '.log', pdf_filename))
        os.unlink(re.sub('\.pdf$', '.aux', pdf_filename))

    def add_text(self, text):
        """Add text to the report."""
        # Write the text to a temporary text file
        text_file = NamedTemporaryFile(suffix='.txt', mode='w')
        text_file.write(text)
        text_file.flush()
        print text_file.name
        self.tempFileList.append(text_file)

    def add_text_file(self, filename):
        """Include the contents of a text file in the report."""
        text_file = open(filename, 'r')
        text = text_file.read()
        self.addText(text)
        text_file.close()

    def add_current_figure(self):
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

    def add_figure(self, figure):
        """Adds the figure to the report.

        Uses the figure method savefig to export the figure to a
        temporary PDF file which is incorporated into the report.
        """

        fig_file = NamedTemporaryFile(suffix='.pdf')
        figure.savefig(fig_file.name)
        #pdf_file = PdfPages(fig_file.name)
        #pdf_file.savefig()
        #pdf_file.close()
        self.tempFileList.append(fig_file)

    def write_report(self, outputfilename='report'):
        """Output the complete report to the given PDF file.
        
        Parameters
        ----------
        outputfilename : string
            The name of the output file for the report. *Note* that the .pdf
            extension will be appended to the given filename.
        """

        output = PdfFileWriter()

        for tempFile in self.tempFileList:
            tempFileName = tempFile.name
            fileType =  tempFileName[(len(tempFileName)-3):len(tempFileName)]
            input = None

            if fileType == 'pdf':
                input = PdfFileReader(file(tempFileName, "rb"))
                os.unlink(tempFileName)
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

