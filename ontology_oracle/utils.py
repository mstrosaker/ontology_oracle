# Copyright (c) 2015 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

import sys, time, urllib2, pandas, math

class DownloadError(Exception):
    def __init__(self, value, tries):
        self.value = value
        self.tries = tries
    def __str__(self):
        s = [repr(self.value)]
        if self.tries > 1:
            s.append('(%s tries)' % self.tries)
        return ' '.join(s)

def download(url, tries=2):
    n_tries = 0
    #print url

    while True:
        n_tries += 1
        try:
            response = urllib2.urlopen(url)
            break
        except urllib2.HTTPError as e:
            if n_tries >= tries:
                raise DownloadError('HTTP response: error %d' % e.code,
                                    n_tries)
            else:
                time.sleep(3)
        except urllib2.URLError as e:
            if n_tries >= tries:
                raise DownloadError('URL error (%d): %s' % (e.reason[0],
                                    e.reason[1]), n_tries)
            else:
                time.sleep(3)
        except:
            if n_tries >= tries:
                raise DownloadError('Unexpected error: %s' % \
                                    sys.exc_info()[0], n_tries)
            else:
                time.sleep(3)

    f = response.read()
    return f

def fold_change(frm, to, log2=True, round_values=0.01):
    """
    Calculates the fold change of "to" as compared to "frm" (from).

    log2: specifies whether the fold change should be log base 2 transformed
    round_values: specifies what value should be added to "fr" and "to" prior
        to calculating the fold change to account for biases
        (see: http://bioinfo.aizeonpublishers.net/content/2013/6/285-292.html)
    """
    frm = frm + round_values
    to = to + round_values

    if frm == 0 and to == 0:
        if log2:
            return 0.0
        return 1.0
    elif frm == 0:
        return float('Inf')
    elif to == 0:
        return -float('Inf')

    if log2:
        return math.log(to, 2) - math.log(frm, 2)
    return to / fr

class dataset:
    def __init__(self, location, format='csv'):
        self.location = location
        self.format = format
        if format == 'csv':
            self.data = pandas.read_csv(location, sep=',')
        elif format.startswith('tab'):
            self.data = pandas.read_table(location)
        self.colnames = self.data.columns.values

    @property
    def rows(self):
        for i, row in self.data.iterrows():
            yield row

    def _make_index(self, idx, frm, to=None):
        for i, row in self.data.iterrows():
            if isinstance(row[frm], basestring) and row[frm] != '':
                if to is None:
                    idx[row[frm]] = row
                else:
                    if isinstance(row[to], basestring):
                        idx[row[frm]] = row[to]

    def index(self, frm, to=None):
        idx = {}
        if isinstance(frm, (list, tuple)):
            for f in frm:
                self._make_index(idx, f, to)
        else:
            self._make_index(idx, frm, to)

        return idx

