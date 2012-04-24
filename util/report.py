from tempfile import NamedTemporaryFile
import txt2pdf
from pyPdf import PdfFileWriter, PdfFileReader
import os
from pylab import savefig
from matplotlib.backends.backend_pdf import PdfPages

# A class to facilitate the creation of model analysis documents
# incorporating both figures (from simulations, etc.) and text.
class Report():
  tempFileList = []

  def addText(self, text):
    # Write the text to a temporary text file
    text_file = NamedTemporaryFile(suffix='.txt', mode='w')
    text_file.write(text)
    text_file.flush()
    print text_file.name
    self.tempFileList.append(text_file)

  def addTextFile(self, filename):
    text_file = open(filename, 'r')
    text = text_file.read()
    addText(text)
    text_file.close()

  def addCurrentFigure(self):
    fig_file = NamedTemporaryFile(suffix='.pdf')
    savefig(fig_file.name)
    #pdf_file = PdfPages(fig_file.name)
    #pdf_file.savefig()
    #pdf_file.close()
    self.tempFileList.append(fig_file)

  def writeReport(self, outputfilename='report'):
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

