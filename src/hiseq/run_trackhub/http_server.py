#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Convert the file_path on server to public http_url
"""


import os
import sys
import re
from hiseq.utils.utils import log, update_obj, run_shell_cmd
from hiseq.utils.file import file_abspath


class HttpServer(object):
    """
    Transform the path of local files to urls (http)

    Attributes:
        s (str): Path to a file/directory
        http_root_dir (str): Path to the http_root,
            it is /var/www/html for apache2
        http_root_url (str): Url to the root of the website,
            default is the IP address

    Example:
    >>> args = {
        's': '/data/public/hub.txt',
        'http_root_dir': '/data/public',
        'http_root_alias': '/public',
        'http_root_url': None,
        'is_https': False
        }
    >>> HttpServer(**args).to_url()
        # for example:
        Alias /public /data/public
        http_root_alias = '/public'
        http_root_dir = '/data/public'
        http_root_url = http://ip/public

    Suppose:
    - Apache2 is properly installed, and working well

    1. Guess the ip address
    2. Get the http_root_dir from argument
    3. Check if the query_path is in the sub_folder of http_root_dir
    4. Transform the query_path to http url

    Parameters
    -----------
    remote_dir : str
        Absolute path to the directory/file in the http_root_dir

    http_root_dir : str
        Absolute path to the root_dir of http config, could be found in file
        `/etc/apache2/sites-available/000-browser.conf`, or other *.conf
        files, if you changed the Apache2 configuration. Check the Alias
        also.

    http_root_url : str
        The URL correspond to the http_root_dir.

    # Default Apache2 configuration

    # Search the http_root_dir:
    On Ubuntu 18.04, could find these config in file:
    `/etc/apache2/sites-available/000-default.conf`, or your *.conf
    $ grep -i 'DocumentRoot' /etc/apache2/sites-available/000-default.conf
        DocumentRoot /var/www/html
    it means, the http_root_dir is "/var/www/html"

    # Fetch the IP of the server
    $ hostname -I
    1.2.3.4 # (fake ip, here)

    # in case multiple hosts, choose first one
    $ hostname -I
    1.2.3.4  5.6.7.8
    you can find the http_root_ip is '1.2.3.4'

    ## For VirtualHost
    If VirtualHost has been build on the server, we can host the files on
    other directories beside the `DocumentRoot /var/www/html`
    for example:

    ```
    # file: /etc/apache2/sites-available/000-default.conf # or other .conf

    Alias /public /data/public
    <Directory "/data/public">
            Require all granted
    </Directory>
    ```

    The `http_root_url` for the virtualHost was: `hostname/publi`, and the
    correspond `http_root_path` was: `/data/public`
    """

    def __init__(
        self,
        s="",
        http_root_dir="",
        http_root_alias="",
        is_https=False,
        **kwargs,
    ):
        self = update_obj(self, kwargs, force=True)
        self.s = s
        self.http_root_dir = http_root_dir
        self.http_root_alias = http_root_alias
        self.is_https = is_https
        self.init_files()

    def init_files(self):
        """
        Initiate the http_root_{dir,alias,url}

        Purpose
        fetch values for: http_root_{dir,alias,url}

        Option-1:
        http_root_dir is a valid URL,
          - parse http_root_alias, string after the {ip|domain},
            http://www.123.com{http_root_alias}

        Option-2:
        http_root_alias if a valid dirname,
          - parse the IP of 'localhost' HTTP server using shell command:
            `hostname -I`
          - determine HTTP server type: 'https' or 'http'
          - construct the http_root_url:
            {http}://{ip}{http_root_alias}
        """
        # absolute path
        self.s = file_abspath(self.s).rstrip("/")
        self.http_root_dir = file_abspath(self.http_root_dir).rstrip("/")

        # update http_root_{alias,url}
        if self.is_url(self.http_root_url):
            self.http_root_url = self.http_root_url.rstrip("/")
            url = re.sub("^((https?)|^(ftp))://", "", self.http_root_url)
            url_ip, url_dir = url.split("/", 1)  #
            self.http_root_alias = "/" + url_dir
        elif len(self.http_root_alias) > 0:
            self.http_root_alias = self.http_root_alias.rstrip("/")
            url_ip = self.get_ip()
            url_http = "https" if self.is_https else "http"
            self.http_root_url = "{http}://{ip}{alias}".format(
                http=url_http, ip=url_ip, alias=self.http_root_alias
            )
        else:
            log.error("check http_root_alias, http_root_url")
            raise ValueError(self.http_root_alias, self.http_root_url)

        # log message
        self.msg = "\n".join(
            [
                "Parameters",
                "{:>16s}: {}, {}".format("s", isinstance(self.s, str), self.s),
                "{:>16s}: {}, {}".format(
                    "http_root_dir",
                    isinstance(self.http_root_dir, str),
                    self.http_root_dir,
                ),
                "{:>16s}: {}, {}".format(
                    "http_root_alias",
                    isinstance(self.http_root_alias, str),
                    self.http_root_alias,
                ),
                "{:>16s}: {}, {}".format(
                    "http_root_url",
                    isinstance(self.http_root_url, str),
                    self.http_root_url,
                ),
            ]
        )

        # check http_root_dir and s
        # s contains http_root_dir
        if not self.s.startswith(self.http_root_dir):
            log.error("s and http_root_dir not match")
            raise ValueError(self.msg)

        # argument type
        if not all(
            [
                isinstance(i, str)
                for i in [
                    self.s,
                    self.http_root_dir,
                    self.http_root_url,
                    self.http_root_alias,
                ]
            ]
        ):
            raise ValueError(self.msg)

    def get_ip(self):
        """
        Extract the public ip address of the server
            - Ubuntu
            $ hostname -I
            1.2.3.4

            $ hostname -I
            1.2.3.4   5.6.7.8
            # in case multiple hosts
        """
        #         return run_shell_cmd('hostname -I')[1].strip()
        ip = run_shell_cmd("hostname -I")[1].strip()
        ip_list = ip.split(" ")
        if len(ip_list) > 1:
            ip = ip_list[0]  # first
            log.warning(
                "Multiple ip addresses found: {}, choose: [{}]".format(
                    ip_list, ip
                )
            )
        # validate
        if self.is_ipv4(ip):
            return ip
        else:
            raise ValueError("not valid ip address found: {}".format(ip))

    def is_ipv4(self, s):
        # validate ipv4
        def is_valid(i):
            try:
                return str(int(i)) == i and 0 <= int(i) <= 255
            except:
                return False

        # check
        return s.count(".") == 3 and all(is_valid(i) for i in s.split("."))

    def is_url(self, s):
        """
        Check if s is url:
        Attributes
            s (str): The url, http://, https://, ftp://
        auto fix: www.abc.com -> http://www.abc.com
        127.0.0.1
        www.name.com
        http://127.0.0.1
        http:name.com
        """
        if not isinstance(s, str):
            return False
        # autofix: add 'http://' to url
        # www.abc.com
        # 127.0.0.1
        p1 = re.compile(
            "(^www)|(^\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})", flags=re.IGNORECASE
        )
        if p1.search(s):
            s = "http://" + s
        # search url
        p2 = re.compile("^((https?)|^ftp)://", flags=re.IGNORECASE)
        return p2.search(s)

    def to_url(self, s=None):
        """
        Convert s -> s_alias -> s_url
        Generate the URL of the hub_txt file, (validate the URL?)
        Attributes
            s (str): Path to the file/directory
        """
        if not isinstance(s, str):
            s = self.s

        if s.startswith(self.http_root_dir):
            s_alias = self.s[len(self.http_root_dir) :]
            return self.http_root_url + s_alias
        else:
            log.error("s and http_root_dir not in the same dir")
            raise ValueError(self.msg)
