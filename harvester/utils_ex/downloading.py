import gzip
import logging
import shutil
import tempfile
import urllib2

import requests
from pathlib import Path


def acquire_files(data_files):
    """
    Ensures that files specified in 'data_files' are the most recent copies.

    The 'data_files' dict should have the following shape:
    {
        <string; filename>: {
            'path': <string; location of the target file on disk>,
            'url': <string; URL of the source (ftp or http/s)>,
            'compressed': <boolean, default False; whether to decompress the file after downloading>
        },
        ...
    }

    :param data_files: a dicts with keys corresponding to the files and values describing how to fetch the files
     (described above)
    :return: A dict with the same keys as 'data_files', but with Path objects as the values
    """

    result = {}

    for name, meta in data_files.items():
        logging.info("Acquiring entry %s from %s..." % (name, meta['url']))

        # unconditionally overwrite each file, since we're downloading it anyway and might as well
        # have the latest copy
        try:
            if meta['url'].startswith('ftp'):
                logging.info("Downloading %(url)s via FTP, storing in %(path)s" % meta)
                r = urllib2.urlopen(meta['url'])
                content = r.read()
            else:
                logging.info("Downloading %(url)s via HTTP, storing in %(path)s" % meta)
                r = requests.get(meta['url'], allow_redirects=True)
                content = r.content

            if meta.has_key('compressed') and meta['compressed']:
                with tempfile.TemporaryFile() as tmp_fp:
                    # use a temporary file to store the compressed content...
                    tmp_fp.write(content)
                    tmp_fp.seek(0)
                    decompressed_fp = gzip.GzipFile(fileobj=tmp_fp, mode='rb')

                    # ...then write it out to its final location while decompressing it
                    with open(meta['path'], 'wb') as f_out:
                        shutil.copyfileobj(decompressed_fp, f_out)
            else:
                # just write the file out directly
                with open(meta['path'], 'wb') as fp:
                    fp.write(content)

            result[name] = Path(meta['path'])

        except requests.exceptions.RequestException as e:
            logging.exception("Couldn't get %s from %s, using cached copy" % (name, meta['url']))

    return result
