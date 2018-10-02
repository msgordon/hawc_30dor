#! /usr/bin/env python
import argparse
from pathlib import Path
from nbconvert import HTMLExporter
import nbformat
from bs4 import BeautifulSoup

HEAD = '''
<!-- Custom stylesheet, it must be in the same directory as the html file -->
<link rel="stylesheet" href="custom.css">
'''

HEAD_REP = '''
<!-- Bootswatch theme -->
<link rel="stylesheet" href="bootstrap.min.css">

<!-- Custom stylesheet, it must be in the same directory as the html file -->
<link rel="stylesheet" href="custom.css">
'''

def trim_bytes(filestring):
    from PIL import Image, ImageChops
    from io import BytesIO
    import base64

    prefix,imstring = filestring.split(',')
    imbytes = base64.b64decode(imstring)
    
    with BytesIO(imbytes) as b:
        image = Image.open(b)
        image.load()

    bg = Image.new(image.mode, image.size, image.getpixel((0,0)))
    diff = ImageChops.difference(image, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()

    new_image = image.crop(bbox)

    with BytesIO() as output:
        new_image.save(output,format='PNG')
        payload = base64.b64encode(output.getvalue()).decode()

    payload = ','.join((prefix,payload))
    return payload


def convert(notebook,outfile=None, title='Notebook',
            hide_warnings=False, hide_info=False,
            trim=False):

    if not outfile:
        outfile = Path(notebook).with_suffix('.html')

    html_exporter = HTMLExporter()
    #with open(notebook, 'r') as f:
    nb = nbformat.read(notebook,as_version=nbformat.NO_CONVERT)

    (body,resources) = html_exporter.from_notebook_node(nb)

    body = body.replace(HEAD,HEAD_REP)
    

    soup = BeautifulSoup(body, 'html.parser')
    soup.head.title.string = title

    if hide_warnings:
        for cell in soup.select('.output_stderr'):
            cell['class'].append('hide')

    if hide_info:
        for cell in soup.select('.output_stdout'):
            if cell.pre and ('INFO: ' in cell.pre.string):
                cell['class'].append('hide')

    if trim:
        for cell in soup.select('.output_area .rendered_html img'):
            cell['src'] = trim_bytes(cell['src'])
    
    with open(outfile, 'wb') as f:
        f.write(soup.encode(formatter=None))

def main():
    parser = argparse.ArgumentParser(description='Convert ipynb to html')
    parser.add_argument('notebook',type=str,help='.ipynb filename')
    parser.add_argument('-o',type=str,metavar='outfile',help='Output filename (default=input with html ext)')
    parser.add_argument('-title',type=str,default='Notebook',help='Meta title for html page (default=Notebook)')
    parser.add_argument('--trim',action='store_true',help='Trim whitespace from figures (requires PIL)')
    parser.add_argument('--hw',action='store_true',help='Hide warning output cells')
    parser.add_argument('--hi',action='store_true',help='Hide logging info cells')
    
    args = parser.parse_args()

    convert(args.notebook, args.o, title=args.title,
            hide_warnings=args.hw, hide_info=args.hi,
            trim=args.trim)



if __name__ == '__main__':
    main()
