"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import

import sys
import os
import click

__author__ = 'Allison MacLeay'


folders_of_interest = [
    os.path.join(os.path.join(os.path.dirname(__file__)), 'tools'),
]
cli_files = {}


msg = '''hi'''


def has_cli_method(script_path):
    """
    Check if a script has a cli() method in order to add it to the main
    :param script_path: to a python script inside Cider packages
    :return: Boolean
    """
    file_obj = open(script_path, 'r').read()
    return "cli()" in file_obj


class MyCLI(click.MultiCommand):
    """This class find .py files that should be part of the CLI provided by click"""
    def list_commands(self, ctx):
        """Parse *.py files and find any script with cli() method"""
        rv_all = []
        for folder in folders_of_interest:
            rv_part = []
            for filename in os.listdir(folder):
                if filename.endswith('.py') and not filename.startswith("__init__"):
                    if not has_cli_method(os.path.join(folder, filename)):
                        continue
                    rv_part.append(filename[:-3])
                    cli_files[filename[:-3]] = folder
            rv_part.sort()
            rv_all.extend(rv_part)  # to sort pipelines then helpers instead of mixing them when help message is printed
        return rv_all

    # pylint: disable=arguments-differ
    def get_command(self, ctx, name):
        """Given a click context returns the command name to be used with the main CLI"""
        ns_all = {}
        if not cli_files:
            self.list_commands(ctx)
        try:
            file_path = os.path.join(cli_files[name], name + '.py')
        except ValueError as err:
            sys.stderr.write('ERROR: Unknown command %s\n' % str(err))
            sys.exit(1)
        with open(file_path) as file_obj:
            # eval is dangerous but the syntax comes from this click library documentation
            # http://click.pocoo.org/5/commands/#custom-multi-commands
            code = compile(file_obj.read(), file_obj.name, 'exec')
            eval(code, ns_all, ns_all)  # pylint: disable=eval-used
        return ns_all['cli']


cli = MyCLI(help=msg)
