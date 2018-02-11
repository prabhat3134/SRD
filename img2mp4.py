#! /usr/bin/python3

import os, sys, re, tempfile, shutil

def main( name, fps=15, *convertargs ):
  regexp = re.compile( name + r'(\d+)[.](jpg|png)$', re.IGNORECASE )
  images = []
  imgext = False if convertargs else None
  for arg in os.listdir( os.curdir ):
    if not arg.startswith( name ):
      continue
    argnum, argext = arg[len(name):].rsplit('.',1)
    if not argnum.isdigit() or argext not in ('jpg','png'):
      print( 'ignoring', arg )
      continue
    if imgext is None:
      imgext = argext
    elif imgext != argext:
      imgext = False
    images.append(( int(argnum), os.path.abspath(arg) ))
  assert images, 'no images found in sequence'
  print( 'found', len(images), 'images' )
  images.sort()
  tmpdir = tempfile.mkdtemp( dir=os.curdir )
  try:
    pattern = os.path.join( tmpdir, 'img%05d.' + (imgext or 'ppm') )
    for i, (n,path) in enumerate( images ):
      if imgext:
        os.symlink( path, pattern % i )
      else:
        print( 'converting', path )
        os.system( 'convert %s %s %s' % ( path, ' '.join(convertargs), pattern % i ) )
    output = name + '.mp4'
    os.system( 'avconv -framerate %s -y -f image2 -i "%s" -r %s -vcodec libx264 -crf 10 %s' % ( fps, pattern, fps, output ) )
  finally:
    print( 'cleaning up..' )
    shutil.rmtree( tmpdir )
  print( 'output written to', output )
  print( 'length: %.2f sec' % ( len(images) / float(fps) ) )
  print( 'file size: %.1f MB' % ( os.path.getsize( output ) / 1024.**2 ) )

if __name__ == '__main__':
  main( *sys.argv[1:] )
