import datetime
import os
from astropy.io import ascii
from astropy.time import Time
import requests
import sys

import logging

SSOS_URL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssos.pl"
RESPONSE_FORMAT = 'tsv'
NEW_LINE = '\r\n'


class ParamDictBuilder(object):
    """ Build a dictionary of parameters needed for an SSOS Query. """

    def __init__(self,
                 mbcobject,
                 verbose=False,
                 search_start_date=Time('2013-01-01', scale='utc'),
                 search_end_date=Time('2017-01-01', scale='utc'),
                 search_method='bynameMPC',
                 error_units='arcseconds',
                 resolve_extension=True,
                 resolve_position=True):
        self.mbcobject = mbcobject
        self.verbose = verbose
        self.search_start_date = search_start_date
        self.search_end_date = search_end_date
        self.search_method = search_method
        self.error_units = error_units
        self.resolve_extension = resolve_extension
        self.resolve_position = resolve_position

    @property
    def mbcobject(self):
        """
        The mbcobjects to be used searched 
        """
        return self._mbcobject

    @mbcobject.setter
    def mbcobject(self, mbcobject):

        self._mbcobject = mbcobject

        """
        self._mbcobject = []
        for mbc in mbcobject:
            assert isinstance(mbcobject, input_mbc_lines)
            if not mbcobject.null_observation:
                self._mbcobject.append(mbc)
        """

    @property
    def verbose(self):
        """
        In verbose mode the SSOS query will return diagnoistic
        information about how the search was done.
        """
        return self._verbose

    @verbose.setter
    def verbose(self, verbose):
        self._verbose = ( verbose and 'yes') or 'no'

    @property
    def search_start_date(self):
        """
        astropy.io.Time object. The start date of SSOS search window.
        """
        return self._search_start_date

    @search_start_date.setter
    def search_start_date(self, search_start_date):
        """
        :type search_start_date: astropy.io.Time
        :param search_start_date: search for frames take after the given date.
        """
        assert isinstance(search_start_date, Time)
        self._search_start_date = search_start_date.replicate(format='iso')
        self._search_start_date.out_subfmt = 'date'

    @property
    def search_end_date(self):
        """
        astropy.io.Time object. The end date of SSOS search window.
        """
        return self._search_end_date

    @search_end_date.setter
    def search_end_date(self, search_end_date):
        """
        :type search_end_date: astropy.io.Time
        :param search_end_date: search for frames take after the given date.
        """
        assert isinstance(search_end_date, Time)
        self._search_end_date = search_end_date.replicate(format='iso')
        self._search_end_date.out_subfmt = 'date'

    @property
    def search_method(self):
        """ 
        When searching by name, select ephemeris
        Must be one of ['bynameCADC', 'bynameMPC', 'bynameLowell', 'bynameHorizons']
        """
        return self._search_method

    @search_method.setter
    def search_method(self, search_method):
        try:
            assert search_method in ['bynameCADC', 'bynameMPC', 'bynameLowell', 'bynameHorizons']
        except AssertionError, e:
            raise AssertionError(
                '{} not in {}'.format(search_method, ['bynameCADC', 'bynameMPC', 'bynameLowell', 'bynameHorizons']))
        self._search_method = search_method

    @property
    def error_units(self):
        """
        For fitting by arc
        """
        return self._error_units

    @error_units.setter
    def error_units(self, error_units):
        """
        :param error_units: 
        """

        self._error_units = error_units

    @property
    def resolve_extension(self):
        """
        Should SSOS resolve and return which extension of a frame the
        object would be in?
        """
        return self._resolve_extension

    @resolve_extension.setter
    def resolve_extension(self, resolve_extension):
        if str(resolve_extension).lower() == "no":
            resolve_extension = False
        self._resolve_extension = (resolve_extension and "yes") or "no"

    @property
    def resolve_position(self):
        """
        Should SSOS resolve and return the predicted X/Y location of
        the source?
        """
        return self._resolve_position

    @resolve_position.setter
    def resolve_position(self, resolve_position):
        if str(resolve_position).lower() == "no":
            resolve_position = False
        self._resolve_position = (resolve_position and "yes") or "no"

    @property
    def params(self):
        """
        The SSOS Query parameters as dictionary, appropriate for url_encoding
        """
        return dict(format=RESPONSE_FORMAT,
                    verbose=self.verbose,
                    epoch1=str(self.search_start_date),
                    epoch2=str(self.search_end_date),
                    search=self.search_method,
                    eunits=self.error_units,
                    extres=self.resolve_extension,
                    xyres=self.resolve_position,
                    object=self.mbcobject
                    )

        # obs=NEW_LINE.join((str(observation) for observation in self.observations))


class Query(object):
    """
    Query the CADC's (or MPC, Lowell, Horizon's) Solar System Object search for a given set of
    main belt comets.

    Inputs:
        - a list of main belt comets

    Optional:
        - a tuple of the start and end times to be searched
          between. Format '%Y-%m-%d'

    Otherwise the temporal range defaults to spanning from the start
    of OSSOS surveying on 2013-01-01 to the present day.

    """

    def __init__(self,
                 mbcobject,
                 search_start_date=Time('2013-01-01', scale='utc'),
                 search_end_date=Time('2017-01-01', scale='utc')):
        self.param_dict_builder = ParamDictBuilder(mbcobject,
                                                   search_start_date=search_start_date,
                                                   search_end_date=search_end_date)

        self.headers = {'User-Agent': 'OSSOS Images'}

    def get(self):
        """
        :return: astropy.table.table
        :raise: AssertionError
        """
        params = self.param_dict_builder.params
        # print("{}\n".format(params))

        self.response = requests.post(SSOS_URL, data=params, headers=self.headers)
        # print(self.response.url)
        try:
            assert isinstance(self.response, requests.Response)
        except AssertionError, e:
            raise AssertionError('response not instance of requests.Response')
        try:
            assert (self.response.status_code == requests.codes.ok )
        except AssertionError, e:
            raise AssertionError('response.status_code =! requests.codes.ok')

        lines = self.response.content
        # note: spelling 'occured' is in SSOIS
        if len(lines) < 2 or "An error occured getting the ephemeris" in lines:
            raise IOError(os.errno.EACCES, "call to SSOIS failed on format error")

        if os.access("backdoor.tsv", os.R_OK):
            lines += open("backdoor.tsv").read()

        return lines