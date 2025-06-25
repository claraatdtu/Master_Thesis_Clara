#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# Author: clsor
# GNU Radio version: 3.10.10.0

from PyQt5 import Qt
from gnuradio import qtgui
from gnuradio import analog
from gnuradio import blocks
from gnuradio import digital
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio import gr
from gnuradio.fft import window
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import soapy
import BPSKRX_wform as wform  # embedded python module
import math
import sip



class BPSKRX(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Not titled yet", catch_exceptions=True)
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Not titled yet")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except BaseException as exc:
            print(f"Qt GUI: Could not set Icon: {str(exc)}", file=sys.stderr)
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "BPSKRX")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Variables
        ##################################################
        self.constellation = constellation = digital.constellation_bpsk().points()
        self.sf_lora = sf_lora = 7
        self.bw_lora = bw_lora = 125000
        self.M = M = len(constellation )
        self.samp_rate = samp_rate = 1000000
        self.bps = bps = int(math.log(M,2))
        self.Rb = Rb = (sf_lora*bw_lora)/2**sf_lora
        self.Sps = Sps = int((bps*samp_rate)/Rb)
        self.ntaps = ntaps = 16*Sps
        self.m = m = int(math.log2(M))
        self.beta = beta = 1
        self.sig_power = sig_power = 1
        self.packet_len = packet_len = 240
        self.num_samples = num_samples = 100000
        self.h = h = wform.rrcos(Sps,ntaps,beta)
        self.delay = delay = bps*8
        self.center_freq = center_freq = 868100000
        self.bw = bw = (2**m)*Rb/m
        self.bpsk = bpsk = digital.constellation_bpsk().base()
        self.bpsk.set_npwr(1.0)
        self.Rs = Rs = Rb/bps
        self.Fmax = Fmax = samp_rate/2

        ##################################################
        # Blocks
        ##################################################

        self.soapy_hackrf_source_0 = None
        dev = 'driver=hackrf'
        stream_args = ''
        tune_args = ['']
        settings = ['']

        self.soapy_hackrf_source_0 = soapy.source(dev, "fc32", 1, 'clock=external',
                                  stream_args, tune_args, settings)
        self.soapy_hackrf_source_0.set_sample_rate(0, 1000000)
        self.soapy_hackrf_source_0.set_bandwidth(0, 0)
        self.soapy_hackrf_source_0.set_frequency(0, center_freq)
        self.soapy_hackrf_source_0.set_gain(0, 'AMP', False)
        self.soapy_hackrf_source_0.set_gain(0, 'LNA', min(max(0, 0.0), 40.0))
        self.soapy_hackrf_source_0.set_gain(0, 'VGA', min(max(16, 0.0), 62.0))
        self.qtgui_time_sink_x_2_0_0_0_1 = qtgui.time_sink_c(
            (32*bps), #size
            samp_rate, #samp_rate
            'AFTER squelch', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0_0_0_1.set_update_time(0.010)
        self.qtgui_time_sink_x_2_0_0_0_1.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_2_0_0_0_1.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_2_0_0_0_1.enable_tags(True)
        self.qtgui_time_sink_x_2_0_0_0_1.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_2_0_0_0_1.enable_autoscale(False)
        self.qtgui_time_sink_x_2_0_0_0_1.enable_grid(False)
        self.qtgui_time_sink_x_2_0_0_0_1.enable_axis_labels(True)
        self.qtgui_time_sink_x_2_0_0_0_1.enable_control_panel(False)
        self.qtgui_time_sink_x_2_0_0_0_1.enable_stem_plot(True)


        labels = ['t4', 'r4', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [3, 3, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [0, 0, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(2):
            if len(labels[i]) == 0:
                if (i % 2 == 0):
                    self.qtgui_time_sink_x_2_0_0_0_1.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_2_0_0_0_1.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_2_0_0_0_1.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_2_0_0_0_1.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_2_0_0_0_1.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_2_0_0_0_1.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_2_0_0_0_1.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_2_0_0_0_1.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_2_0_0_0_1_win = sip.wrapinstance(self.qtgui_time_sink_x_2_0_0_0_1.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_2_0_0_0_1_win)
        self.qtgui_time_sink_x_2_0_0_0_0_0_0 = qtgui.time_sink_c(
            (32*bps), #size
            samp_rate, #samp_rate
            'after agc', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_update_time(0.010)
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_2_0_0_0_0_0_0.enable_tags(True)
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.enable_grid(False)
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.enable_stem_plot(True)


        labels = ['t4', 'r4', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [3, 3, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [0, 0, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(2):
            if len(labels[i]) == 0:
                if (i % 2 == 0):
                    self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_2_0_0_0_0_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_2_0_0_0_0_0_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_2_0_0_0_0_0_0_win)
        self.qtgui_time_sink_x_2_0_0_0_0_0 = qtgui.time_sink_c(
            (32*bps), #size
            samp_rate, #samp_rate
            'after clock reco', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0_0_0_0_0.set_update_time(0.010)
        self.qtgui_time_sink_x_2_0_0_0_0_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_2_0_0_0_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_2_0_0_0_0_0.enable_tags(True)
        self.qtgui_time_sink_x_2_0_0_0_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_2_0_0_0_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_2_0_0_0_0_0.enable_grid(False)
        self.qtgui_time_sink_x_2_0_0_0_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_2_0_0_0_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_2_0_0_0_0_0.enable_stem_plot(True)


        labels = ['t4', 'r4', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [3, 3, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [0, 0, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(2):
            if len(labels[i]) == 0:
                if (i % 2 == 0):
                    self.qtgui_time_sink_x_2_0_0_0_0_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_2_0_0_0_0_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_2_0_0_0_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_2_0_0_0_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_2_0_0_0_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_2_0_0_0_0_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_2_0_0_0_0_0_win)
        self.qtgui_time_sink_x_2_0_0_0_0 = qtgui.time_sink_c(
            (32*bps), #size
            samp_rate, #samp_rate
            'after const', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0_0_0_0.set_update_time(0.010)
        self.qtgui_time_sink_x_2_0_0_0_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_2_0_0_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_2_0_0_0_0.enable_tags(True)
        self.qtgui_time_sink_x_2_0_0_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_2_0_0_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_2_0_0_0_0.enable_grid(False)
        self.qtgui_time_sink_x_2_0_0_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_2_0_0_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_2_0_0_0_0.enable_stem_plot(True)


        labels = ['t4', 'r4', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [3, 3, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [0, 0, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(2):
            if len(labels[i]) == 0:
                if (i % 2 == 0):
                    self.qtgui_time_sink_x_2_0_0_0_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_2_0_0_0_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_2_0_0_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_2_0_0_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_2_0_0_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_2_0_0_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_2_0_0_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_2_0_0_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_2_0_0_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_2_0_0_0_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_2_0_0_0_0_win)
        self.qtgui_time_sink_x_2_0_0_0 = qtgui.time_sink_c(
            (32*bps), #size
            samp_rate, #samp_rate
            'BEFORE squelch', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0_0_0.set_update_time(0.010)
        self.qtgui_time_sink_x_2_0_0_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_2_0_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_2_0_0_0.enable_tags(True)
        self.qtgui_time_sink_x_2_0_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_2_0_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_2_0_0_0.enable_grid(False)
        self.qtgui_time_sink_x_2_0_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_2_0_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_2_0_0_0.enable_stem_plot(True)


        labels = ['t4', 'r4', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [3, 3, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [0, 0, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(2):
            if len(labels[i]) == 0:
                if (i % 2 == 0):
                    self.qtgui_time_sink_x_2_0_0_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_2_0_0_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_2_0_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_2_0_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_2_0_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_2_0_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_2_0_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_2_0_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_2_0_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_2_0_0_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_2_0_0_0_win)
        self.qtgui_time_sink_x_2_0 = qtgui.time_sink_f(
            (32*bps), #size
            samp_rate, #samp_rate
            'Comparison bits', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0.set_update_time(0.010)
        self.qtgui_time_sink_x_2_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_2_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_2_0.enable_tags(True)
        self.qtgui_time_sink_x_2_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_2_0.enable_autoscale(False)
        self.qtgui_time_sink_x_2_0.enable_grid(False)
        self.qtgui_time_sink_x_2_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_2_0.enable_control_panel(False)
        self.qtgui_time_sink_x_2_0.enable_stem_plot(True)


        labels = ['t4', 'r4', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [3, 3, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [0, 0, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(1):
            if len(labels[i]) == 0:
                self.qtgui_time_sink_x_2_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_2_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_2_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_2_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_2_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_2_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_2_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_2_0_win = sip.wrapinstance(self.qtgui_time_sink_x_2_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_2_0_win)
        self.qtgui_sink_x_0_0_0 = qtgui.sink_f(
            1024, #fftsize
            window.WIN_BLACKMAN_hARRIS, #wintype
            center_freq, #fc
            (bw*8), #bw
            'constell', #name
            True, #plotfreq
            True, #plotwaterfall
            True, #plottime
            True, #plotconst
            None # parent
        )
        self.qtgui_sink_x_0_0_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_0_0_win = sip.wrapinstance(self.qtgui_sink_x_0_0_0.qwidget(), Qt.QWidget)

        self.qtgui_sink_x_0_0_0.enable_rf_freq(True)

        self.top_layout.addWidget(self._qtgui_sink_x_0_0_0_win)
        self.interp_fir_filter_xxx_0_0 = filter.interp_fir_filter_ccf(1, h)
        self.interp_fir_filter_xxx_0_0.declare_sample_delay(0)
        self.digital_probe_mpsk_snr_est_c_0 = digital.probe_mpsk_snr_est_c(0, 10000, 0.001)
        self.digital_constellation_decoder_cb_0 = digital.constellation_decoder_cb(bpsk)
        self.digital_clock_recovery_mm_xx_0_0 = digital.clock_recovery_mm_cc((Sps*(1+0.0)), 0.01, 0.5, 0.1, 0.05)
        self.blocks_unpacked_to_packed_xx_0 = blocks.unpacked_to_packed_bb(bps, gr.GR_MSB_FIRST)
        self.blocks_unpack_k_bits_bb_0 = blocks.unpack_k_bits_bb(8)
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_cc(1)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_char*1, 'C:\\Users\\clsor\\OneDrive\\Documents\\MATLAB\\Master_Thesis_Clara\\Master_Thesis_Clara\\3-GNU radio implementation\\SDR files of bits\\QPSKsdroutput_bytes', False)
        self.blocks_file_sink_0.set_unbuffered(False)
        self.blocks_char_to_float_0_1_1 = blocks.char_to_float(1, 1)
        self.blocks_char_to_float_0_1 = blocks.char_to_float(1, 1)
        self.analog_pwr_squelch_xx_0 = analog.pwr_squelch_cc((-20), (1e-4), 0, True)
        self.analog_agc3_xx_0 = analog.agc3_cc((1e-3), (1e-4), 1, 0.5, 1, 1)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_agc3_xx_0, 0), (self.digital_constellation_decoder_cb_0, 0))
        self.connect((self.analog_agc3_xx_0, 0), (self.qtgui_time_sink_x_2_0_0_0_0_0_0, 0))
        self.connect((self.analog_pwr_squelch_xx_0, 0), (self.interp_fir_filter_xxx_0_0, 0))
        self.connect((self.analog_pwr_squelch_xx_0, 0), (self.qtgui_time_sink_x_2_0_0_0_1, 0))
        self.connect((self.blocks_char_to_float_0_1, 0), (self.qtgui_time_sink_x_2_0, 0))
        self.connect((self.blocks_char_to_float_0_1_1, 0), (self.qtgui_sink_x_0_0_0, 0))
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.digital_clock_recovery_mm_xx_0_0, 0))
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.qtgui_time_sink_x_2_0_0_0_0, 0))
        self.connect((self.blocks_unpack_k_bits_bb_0, 0), (self.blocks_char_to_float_0_1, 0))
        self.connect((self.blocks_unpack_k_bits_bb_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.blocks_unpacked_to_packed_xx_0, 0), (self.blocks_unpack_k_bits_bb_0, 0))
        self.connect((self.digital_clock_recovery_mm_xx_0_0, 0), (self.analog_agc3_xx_0, 0))
        self.connect((self.digital_clock_recovery_mm_xx_0_0, 0), (self.qtgui_time_sink_x_2_0_0_0_0_0, 0))
        self.connect((self.digital_constellation_decoder_cb_0, 0), (self.blocks_char_to_float_0_1_1, 0))
        self.connect((self.digital_constellation_decoder_cb_0, 0), (self.blocks_unpacked_to_packed_xx_0, 0))
        self.connect((self.interp_fir_filter_xxx_0_0, 0), (self.blocks_multiply_const_vxx_0, 0))
        self.connect((self.soapy_hackrf_source_0, 0), (self.analog_pwr_squelch_xx_0, 0))
        self.connect((self.soapy_hackrf_source_0, 0), (self.digital_probe_mpsk_snr_est_c_0, 0))
        self.connect((self.soapy_hackrf_source_0, 0), (self.qtgui_time_sink_x_2_0_0_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "BPSKRX")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_constellation(self):
        return self.constellation

    def set_constellation(self, constellation):
        self.constellation = constellation
        self.set_M(len(self.constellation ))

    def get_sf_lora(self):
        return self.sf_lora

    def set_sf_lora(self, sf_lora):
        self.sf_lora = sf_lora
        self.set_Rb((self.sf_lora*self.bw_lora)/2**self.sf_lora)

    def get_bw_lora(self):
        return self.bw_lora

    def set_bw_lora(self, bw_lora):
        self.bw_lora = bw_lora
        self.set_Rb((self.sf_lora*self.bw_lora)/2**self.sf_lora)

    def get_M(self):
        return self.M

    def set_M(self, M):
        self.M = M
        self.set_bps(int(math.log(self.M,2)))
        self.set_m(int(math.log2(self.M)))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_Fmax(self.samp_rate/2)
        self.set_Sps(int((self.bps*self.samp_rate)/self.Rb))
        self.qtgui_time_sink_x_2_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_2_0_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_2_0_0_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_2_0_0_0_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_2_0_0_0_0_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_2_0_0_0_1.set_samp_rate(self.samp_rate)

    def get_bps(self):
        return self.bps

    def set_bps(self, bps):
        self.bps = bps
        self.set_Rs(self.Rb/self.bps)
        self.set_Sps(int((self.bps*self.samp_rate)/self.Rb))
        self.set_delay(self.bps*8)

    def get_Rb(self):
        return self.Rb

    def set_Rb(self, Rb):
        self.Rb = Rb
        self.set_Rs(self.Rb/self.bps)
        self.set_Sps(int((self.bps*self.samp_rate)/self.Rb))
        self.set_bw((2**self.m)*self.Rb/self.m)

    def get_Sps(self):
        return self.Sps

    def set_Sps(self, Sps):
        self.Sps = Sps
        self.set_h(wform.rrcos(self.Sps,self.ntaps,self.beta))
        self.set_ntaps(16*self.Sps)
        self.digital_clock_recovery_mm_xx_0_0.set_omega((self.Sps*(1+0.0)))

    def get_ntaps(self):
        return self.ntaps

    def set_ntaps(self, ntaps):
        self.ntaps = ntaps
        self.set_h(wform.rrcos(self.Sps,self.ntaps,self.beta))

    def get_m(self):
        return self.m

    def set_m(self, m):
        self.m = m
        self.set_bw((2**self.m)*self.Rb/self.m)

    def get_beta(self):
        return self.beta

    def set_beta(self, beta):
        self.beta = beta
        self.set_h(wform.rrcos(self.Sps,self.ntaps,self.beta))

    def get_sig_power(self):
        return self.sig_power

    def set_sig_power(self, sig_power):
        self.sig_power = sig_power

    def get_packet_len(self):
        return self.packet_len

    def set_packet_len(self, packet_len):
        self.packet_len = packet_len

    def get_num_samples(self):
        return self.num_samples

    def set_num_samples(self, num_samples):
        self.num_samples = num_samples

    def get_h(self):
        return self.h

    def set_h(self, h):
        self.h = h
        self.interp_fir_filter_xxx_0_0.set_taps(self.h)

    def get_delay(self):
        return self.delay

    def set_delay(self, delay):
        self.delay = delay

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.qtgui_sink_x_0_0_0.set_frequency_range(self.center_freq, (self.bw*8))
        self.soapy_hackrf_source_0.set_frequency(0, self.center_freq)

    def get_bw(self):
        return self.bw

    def set_bw(self, bw):
        self.bw = bw
        self.qtgui_sink_x_0_0_0.set_frequency_range(self.center_freq, (self.bw*8))

    def get_bpsk(self):
        return self.bpsk

    def set_bpsk(self, bpsk):
        self.bpsk = bpsk
        self.digital_constellation_decoder_cb_0.set_constellation(self.bpsk)

    def get_Rs(self):
        return self.Rs

    def set_Rs(self, Rs):
        self.Rs = Rs

    def get_Fmax(self):
        return self.Fmax

    def set_Fmax(self, Fmax):
        self.Fmax = Fmax




def main(top_block_cls=BPSKRX, options=None):

    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()

    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    qapp.exec_()

if __name__ == '__main__':
    main()
