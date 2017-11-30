#!/bin/bash

set -e


#ignored-modules=cosmos since pylint keep compalining about not being able to import cosmos
pylint --disable=I0011,R0201,E711,E712,E0012 celltics/


# Disable F401: unused imports. Only really applies to __init__.py files, since all other unused imports
# Disable W391:  blank line at end of file but pylint requires a blank line at EOF!

flake8 --max-line-length=120 --ignore=F401,E127,W391,E711,E712 celltics/
