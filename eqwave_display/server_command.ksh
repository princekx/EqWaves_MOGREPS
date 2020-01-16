#!/usr/bin/ksh

netstat -tupln

module load scitools

/opt/scitools/environments/default/2019_02_27/bin/bokeh serve --port 8080 main.py
