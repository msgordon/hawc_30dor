#! /usr/bin/env python
import argparse
from pathlib import Path
from nbconvert import HTMLExporter
import nbformat
from bs4 import BeautifulSoup
from copy import copy
from urllib.request import urlopen

HEAD = '''
<!-- Custom stylesheet, it must be in the same directory as the html file -->
<link rel="stylesheet" href="custom.css">
'''

HEAD_REP = '''
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="shortcut icon" href="https://www.sofia.usra.edu/sites/default/files/favicon.ico" type="image/vnd.microsoft.icon" />

<!-- Bootswatch theme -->
<link rel="stylesheet" href="bootstrap.min.css">

<!-- Custom stylesheet, it must be in the same directory as the html file -->
<link rel="stylesheet" href="custom.css">
'''

HEAD_REM = '''
div.output_area .rendered_html img {
  margin-left: 0;
  margin-right: 0;
}
div.output_area img,
div.output_area svg {
  max-width: 100%;
  height: auto;
}
'''

NAV = '''
<nav class="navbar navbar-inverse navbar-sofia">
  <div class="container">
    <!-- Brand and toggle get grouped for better mobile display -->
    <div class="navbar-header" id="site-title-head">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar-collapse-1" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
      <a class="navbar-brand" href="https://www.sofia.usra.edu/">
        <img alt="SOFIA Science Center" src="https://www.sofia.usra.edu/sites/default/files/sofialogo-blkbg-trans.png">
      </a>
    </div>

    <!-- Collect the nav links, forms, and other content for toggling -->
    <div class="collapse navbar-collapse" id="navbar-collapse-1">
      <ul class="nav navbar-nav navbar-left">
        <li class="site-title no-hover" id="site-li-head">
	  <a href="https://www.sofia.usra.edu" class="no-hover">SOFIA Science Center</a>
	  <a id="site-description" href="https://www.sofia.usra.edu" class="no-hover">
	    Stratospheric Observatory for Infrared Astronomy
	  </a>
	</li>
      </ul>
      <ul class="nav navbar-nav navbar-right">
        <li><a href="https://www.sofia.usra.edu">Home</a></li>
	<li><a href="https://www.sofia.usra.edu/science">For Researchers</a></li>
	<li><a href="https://www.sofia.usra.edu/multimedia">Multimedia</a></li>
      </ul>
    </div><!-- /.navbar-collapse -->
  </div><!-- /.container -->
</nav>
'''

SUBNAV = '''
  <nav class="navbar navbar-inverse navbar-sofia hidden-xs hidden-ms">
  <div class="container-fluid">
    <!-- Brand and toggle get grouped for better mobile display -->
    <div class="navbar-header">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar-collapse-2" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
    </div>

    <!-- Collect the nav links, forms, and other content for toggling -->
    <div class="collapse navbar-collapse" id="navbar-collapse-2">
      <ul class="nav navbar-nav" id="nice-menu">
      {% subnav %}
      </ul>
    </div><!-- /.navbar-collapse -->
  </div><!-- /.container -->
</nav>
'''

JSRC = '''
<!-- jQuery (necessary for Bootstrap's JavaScript plugins) -->
						    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>
						    <!-- Latest compiled and minified JavaScript -->
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>
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
            trim=False, headers=False, copycss=False):

    if not outfile:
        outfile = Path(notebook).with_suffix('.html')

    html_exporter = HTMLExporter()
    
    nb = nbformat.read(notebook,as_version=nbformat.NO_CONVERT)

    (body,resources) = html_exporter.from_notebook_node(nb)

    body = body.replace(HEAD,HEAD_REP)
    body = body.replace(HEAD_REM,'')

    if copycss:
        with open('bootstrap.min.css') as f:
            style = '<style type="text/css">\n%s\n</style>' % f.read()
            body = body.replace('<link rel="stylesheet" href="bootstrap.min.css">',
                                style)
        with open('custom.css') as f:
            style = '<style type="text/css">\n%s\n</style>' % f.read()
            body = body.replace('<link rel="stylesheet" href="custom.css">',
                                style)

        

    soup = BeautifulSoup(body, 'html.parser')
    soup.head.title.string = title

    soup.find(id='notebook').insert_before(NAV)
    soup.body.append(JSRC)

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

    # make img responsive
    for cell in soup.select('.output_area .rendered_html img'):
        cell['class'] = ['img-responsive', 'center-block']
        #del cell['width']

        # make this invisible on xs devices
        divp = cell.find_parent('div')
        divp['class'].append('hidden-xs')
        divp['class'].append('visible-ms')

        #prompt
        prompt = divp.find_previous('div')
        prompt['class'].append('hidden-xs')
        prompt['class'].append('visible-ms')
        
        # add copy that only displays on xs devices
        divp.insert_after(copy(cell))
        newcell = cell.find_next('img')
        newcell['class'] = ['img-responsive', 'center-block','visible-xs-block','hidden-ms']
        del newcell['width']

    if headers:
        cells = get_headers('https://www.sofia.usra.edu/science')

        cells = [str(tag) for tag in cells]
        subnav = SUBNAV.replace('{% subnav %}','\n'.join(cells))

        soup.select('#notebook-container')[0].div.insert_before(subnav)

        '''
        # must write out and read back to add header?
        with open(outfile, 'wb') as f:
            f.write(soup.encode(formatter=None))

        with open(outfile,'r') as f:
            soup = BeautifulSoup(f, 'html.parser')

        print(cells)
        print(soup.select('#nice-menu'))
        soup.select('#nice-menu')[0].append(cells)
        print(soup.select('#nice-menu'))
        '''
    
    with open(outfile, 'wb') as f:
        f.write(soup.encode(formatter=None))


def get_headers(url,base='https://www.sofia.usra.edu'):
    with urlopen(url) as page:
        soup = BeautifulSoup(page,'html.parser')

    cells = []
    for cell in soup.select('ul.nice-menu-down li.menuparent'):
        cell['class'] = 'dropdown'
        cell.a['href'] = '%s%s' % (base, cell.a['href'])
        cell.a['class'] = 'dropdown-toggle'
        cell.a['data-toggle'] = 'dropdown'
        cell.a['role'] = 'button'
        cell.a['aria-haspopup'] = 'true'
        cell.a['aria-expanded'] = 'false'
        cell.ul['class'] = 'dropdown-menu'
        for li in cell.ul.find_all('li'):
            del li['class']
            li.a['href'] = '%s%s' % (base, li.a['href'])

        cells.append(cell)
    return cells

def main():
    parser = argparse.ArgumentParser(description='Convert ipynb to html')
    parser.add_argument('notebook',type=str,help='.ipynb filename')
    parser.add_argument('-o',type=str,metavar='outfile',help='Output filename (default=input with html ext)')
    parser.add_argument('-title',type=str,default='Notebook',help='Meta title for html page (default=Notebook)')
    parser.add_argument('--trim',action='store_true',help='Trim whitespace from figures (requires PIL)')
    parser.add_argument('--hw',action='store_true',help='Hide warning output cells')
    parser.add_argument('--hi',action='store_true',help='Hide logging info cells')
    parser.add_argument('--headers',action='store_true',help='Add header bar from sofia website')
    parser.add_argument('--copycss',action='store_true',help='Copy CSS into header')
    
    args = parser.parse_args()

    convert(args.notebook, args.o, title=args.title,
            hide_warnings=args.hw, hide_info=args.hi,
            trim=args.trim,headers=args.headers,
            copycss=args.copycss)



if __name__ == '__main__':
    main()
