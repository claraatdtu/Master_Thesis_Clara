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
from PyQt5 import QtCore
from gnuradio import analog
import math
from gnuradio import blocks
import pmt
from gnuradio import digital
from gnuradio import fec
from gnuradio import gr
from gnuradio.filter import firdes
from gnuradio.fft import window
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
import BFSKloopbacktry_epy_block_0_0 as epy_block_0_0  # embedded python block
import sip



class BFSKloopbacktry(gr.top_block, Qt.QWidget):

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

        self.settings = Qt.QSettings("GNU Radio", "BFSKloopbacktry")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Variables
        ##################################################
        self.sf_lora = sf_lora = 7
        self.bw_lora = bw_lora = 125000
        self.M = M = 2
        self.m = m = int(math.log2(M))
        self.Rb = Rb = (sf_lora*bw_lora)/2**sf_lora
        self.samp_rate = samp_rate = 1000000
        self.eb_n0_dB = eb_n0_dB = -5
        self.bw = bw = (2**m)*Rb/m
        self.bps = bps = int(math.log(M,2))
        self.sig_power = sig_power = 1
        self.packet_len = packet_len = 240
        self.num_samples = num_samples = 100000
        self.noise_power = noise_power = 1 / (2 * m* 10**(eb_n0_dB / 10))
        self.ndisp = ndisp = 2000
        self.fsk_deviation = fsk_deviation = bw/m
        self.center_freq = center_freq = 868100000
        self.Sps = Sps = int((bps*samp_rate)/Rb)
        self.Rs = Rs = Rb/bps

        ##################################################
        # Blocks
        ##################################################

        self.qtgui_time_sink_x_2_0_0_2_0 = qtgui.time_sink_f(
            (32*bps), #size
            samp_rate, #samp_rate
            'TX float to char 0 1 ', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0_0_2_0.set_update_time(0.10)
        self.qtgui_time_sink_x_2_0_0_2_0.set_y_axis(-0, 3)

        self.qtgui_time_sink_x_2_0_0_2_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_2_0_0_2_0.enable_tags(True)
        self.qtgui_time_sink_x_2_0_0_2_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_2_0_0_2_0.enable_autoscale(False)
        self.qtgui_time_sink_x_2_0_0_2_0.enable_grid(False)
        self.qtgui_time_sink_x_2_0_0_2_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_2_0_0_2_0.enable_control_panel(False)
        self.qtgui_time_sink_x_2_0_0_2_0.enable_stem_plot(True)


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
                self.qtgui_time_sink_x_2_0_0_2_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_2_0_0_2_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_2_0_0_2_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_2_0_0_2_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_2_0_0_2_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_2_0_0_2_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_2_0_0_2_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_2_0_0_2_0_win = sip.wrapinstance(self.qtgui_time_sink_x_2_0_0_2_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_2_0_0_2_0_win)
        self.qtgui_time_sink_x_2_0_0 = qtgui.time_sink_f(
            (32*bps), #size
            samp_rate, #samp_rate
            'Compa bit', #name
            2, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0_0.set_update_time(0.10)
        self.qtgui_time_sink_x_2_0_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_2_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_2_0_0.enable_tags(True)
        self.qtgui_time_sink_x_2_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_2_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_2_0_0.enable_grid(False)
        self.qtgui_time_sink_x_2_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_2_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_2_0_0.enable_stem_plot(True)


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
                self.qtgui_time_sink_x_2_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_2_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_2_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_2_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_2_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_2_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_2_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_2_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_2_0_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_2_0_0_win)
        self.qtgui_time_sink_x_0_2_0 = qtgui.time_sink_f(
            ndisp, #size
            samp_rate, #samp_rate
            'Tx Bits', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_0_2_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0_2_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_0_2_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0_2_0.enable_tags(True)
        self.qtgui_time_sink_x_0_2_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0_2_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0_2_0.enable_grid(True)
        self.qtgui_time_sink_x_0_2_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0_2_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0_2_0.enable_stem_plot(False)

        self.qtgui_time_sink_x_0_2_0.disable_legend()

        labels = ['Signal 1', 'Signal 2', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(1):
            if len(labels[i]) == 0:
                self.qtgui_time_sink_x_0_2_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_0_2_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_2_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_2_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_2_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_2_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_2_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_2_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0_2_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_2_0_win)
        self.qtgui_time_sink_x_0_0 = qtgui.time_sink_f(
            ndisp, #size
            samp_rate, #samp_rate
            'Rx Bits', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_0_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0_0.enable_tags(True)
        self.qtgui_time_sink_x_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0_0.enable_grid(True)
        self.qtgui_time_sink_x_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0_0.enable_stem_plot(False)

        self.qtgui_time_sink_x_0_0.disable_legend()

        labels = ['Signal 1', 'Signal 2', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(1):
            if len(labels[i]) == 0:
                self.qtgui_time_sink_x_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_0_win)
        self.qtgui_time_sink_x_0 = qtgui.time_sink_f(
            1024, #size
            samp_rate, #samp_rate
            'after noise', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0.set_y_axis(-2, 2)

        self.qtgui_time_sink_x_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0.enable_tags(True)
        self.qtgui_time_sink_x_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0.enable_grid(False)
        self.qtgui_time_sink_x_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0.enable_stem_plot(False)


        labels = ['Signal 1', 'Signal 2', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(1):
            if len(labels[i]) == 0:
                self.qtgui_time_sink_x_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_win)
        self.qtgui_sink_x_0_0_0 = qtgui.sink_c(
            1024, #fftsize
            window.WIN_BLACKMAN_hARRIS, #wintype
            center_freq, #fc
            (bw*8), #bw
            'before noise', #name
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
        self.qtgui_number_sink_0_0 = qtgui.number_sink(
            gr.sizeof_float,
            0,
            qtgui.NUM_GRAPH_HORIZ,
            1,
            None # parent
        )
        self.qtgui_number_sink_0_0.set_update_time(0.10)
        self.qtgui_number_sink_0_0.set_title('BER between Tx and Rx')

        labels = ['', '', '', '', '',
            '', '', '', '', '']
        units = ['', '', '', '', '',
            '', '', '', '', '']
        colors = [("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"),
            ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black")]
        factor = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]

        for i in range(1):
            self.qtgui_number_sink_0_0.set_min(i, -1)
            self.qtgui_number_sink_0_0.set_max(i, 1)
            self.qtgui_number_sink_0_0.set_color(i, colors[i][0], colors[i][1])
            if len(labels[i]) == 0:
                self.qtgui_number_sink_0_0.set_label(i, "Data {0}".format(i))
            else:
                self.qtgui_number_sink_0_0.set_label(i, labels[i])
            self.qtgui_number_sink_0_0.set_unit(i, units[i])
            self.qtgui_number_sink_0_0.set_factor(i, factor[i])

        self.qtgui_number_sink_0_0.enable_autoscale(False)
        self._qtgui_number_sink_0_0_win = sip.wrapinstance(self.qtgui_number_sink_0_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_number_sink_0_0_win)
        self.fec_ber_bf_0_0 = fec.ber_bf(False, 100, -7.0)
        self._eb_n0_dB_range = qtgui.Range(-5, 15, 1/100, -5, 200)
        self._eb_n0_dB_win = qtgui.RangeWidget(self._eb_n0_dB_range, self.set_eb_n0_dB, "Eb/N0", "counter_slider", float, QtCore.Qt.Horizontal)
        self.top_grid_layout.addWidget(self._eb_n0_dB_win, 2, 0, 1, 1)
        for r in range(2, 3):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.digital_chunks_to_symbols_xx_0 = digital.chunks_to_symbols_bf([0, 1], 1)
        self.blocks_unpack_k_bits_bb_1 = blocks.unpack_k_bits_bb(m)
        self.blocks_throttle2_0_0 = blocks.throttle( gr.sizeof_gr_complex*1, samp_rate, True, 0 if "auto" == "auto" else max( int(float(0.1) * samp_rate) if "auto" == "time" else int(0.1), 1) )
        self.blocks_stream_to_tagged_stream_0 = blocks.stream_to_tagged_stream(gr.sizeof_char, 1, packet_len, "packet_len")
        self.blocks_skiphead_0_1 = blocks.skiphead(gr.sizeof_char*1, 0)
        self.blocks_skiphead_0_0_0 = blocks.skiphead(gr.sizeof_char*1, 0)
        self.blocks_repeat_1 = blocks.repeat(gr.sizeof_float*1, Sps)
        self.blocks_repeat_0_0_0_0 = blocks.repeat(gr.sizeof_char*1, Sps)
        self.blocks_repeat_0_0_0 = blocks.repeat(gr.sizeof_char*1, Sps)
        self.blocks_pack_k_bits_bb_0 = blocks.pack_k_bits_bb(m)
        self.blocks_multiply_const_vxx_0_0 = blocks.multiply_const_cc(samp_rate/(2*math.pi*fsk_deviation))
        self.blocks_keep_one_in_n_0 = blocks.keep_one_in_n(gr.sizeof_float*1, Sps)
        self.blocks_float_to_uchar_0 = blocks.float_to_uchar(1, 1, 0)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_char*1, 'C:\\Users\\clsor\\OneDrive\\Documents\\MATLAB\\Master_Thesis_Clara\\Master_Thesis_Clara\\3-GNU radio implementation\\SDR files of bits\\sdrinput', False, 0, 0)
        self.blocks_file_source_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_char*1, 'C:\\Users\\clsor\\OneDrive\\Documents\\MATLAB\\Master_Thesis_Clara\\Master_Thesis_Clara\\3-GNU radio implementation\\SDR files of bits\\BFSKsdroutput-5', False)
        self.blocks_file_sink_0.set_unbuffered(False)
        self.blocks_delay_1_0_0_0_0 = blocks.delay(gr.sizeof_char*1, 0)
        self.blocks_char_to_float_0_2_0_0 = blocks.char_to_float(1, 1)
        self.blocks_char_to_float_0_2_0 = blocks.char_to_float(1, 1)
        self.blocks_char_to_float_0_1_0 = blocks.char_to_float(1, 1)
        self.blocks_char_to_float_0_0_1_0 = blocks.char_to_float(1, 1)
        self.blocks_add_xx_1_0 = blocks.add_vcc(1)
        self.analog_quadrature_demod_cf_0_0 = analog.quadrature_demod_cf((samp_rate/(2*math.pi*fsk_deviation)))
        self.analog_noise_source_x_0_0 = analog.noise_source_c(analog.GR_GAUSSIAN, math.sqrt(noise_power), 0)
        self.analog_frequency_modulator_fc_0 = analog.frequency_modulator_fc(((2*math.pi*fsk_deviation)/samp_rate))


        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_frequency_modulator_fc_0, 0), (self.blocks_multiply_const_vxx_0_0, 0))
        self.connect((self.analog_noise_source_x_0_0, 0), (self.blocks_add_xx_1_0, 0))
        self.connect((self.analog_quadrature_demod_cf_0_0, 0), (self.blocks_keep_one_in_n_0, 0))
        self.connect((self.analog_quadrature_demod_cf_0_0, 0), (self.qtgui_time_sink_x_0, 0))
        self.connect((self.blocks_add_xx_1_0, 0), (self.analog_quadrature_demod_cf_0_0, 0))
        self.connect((self.blocks_char_to_float_0_0_1_0, 0), (self.qtgui_time_sink_x_2_0_0, 1))
        self.connect((self.blocks_char_to_float_0_1_0, 0), (self.qtgui_time_sink_x_2_0_0, 0))
        self.connect((self.blocks_char_to_float_0_2_0, 0), (self.qtgui_time_sink_x_0_0, 0))
        self.connect((self.blocks_char_to_float_0_2_0_0, 0), (self.qtgui_time_sink_x_0_2_0, 0))
        self.connect((self.blocks_delay_1_0_0_0_0, 0), (self.blocks_skiphead_0_1, 0))
        self.connect((self.blocks_file_source_0, 0), (self.blocks_delay_1_0_0_0_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.blocks_repeat_0_0_0_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.blocks_stream_to_tagged_stream_0, 0))
        self.connect((self.blocks_float_to_uchar_0, 0), (self.blocks_unpack_k_bits_bb_1, 0))
        self.connect((self.blocks_keep_one_in_n_0, 0), (self.blocks_float_to_uchar_0, 0))
        self.connect((self.blocks_multiply_const_vxx_0_0, 0), (self.blocks_throttle2_0_0, 0))
        self.connect((self.blocks_pack_k_bits_bb_0, 0), (self.digital_chunks_to_symbols_xx_0, 0))
        self.connect((self.blocks_repeat_0_0_0, 0), (self.blocks_char_to_float_0_2_0, 0))
        self.connect((self.blocks_repeat_0_0_0_0, 0), (self.blocks_char_to_float_0_2_0_0, 0))
        self.connect((self.blocks_repeat_1, 0), (self.analog_frequency_modulator_fc_0, 0))
        self.connect((self.blocks_skiphead_0_0_0, 0), (self.blocks_char_to_float_0_0_1_0, 0))
        self.connect((self.blocks_skiphead_0_0_0, 0), (self.fec_ber_bf_0_0, 1))
        self.connect((self.blocks_skiphead_0_1, 0), (self.blocks_char_to_float_0_1_0, 0))
        self.connect((self.blocks_skiphead_0_1, 0), (self.fec_ber_bf_0_0, 0))
        self.connect((self.blocks_stream_to_tagged_stream_0, 0), (self.blocks_pack_k_bits_bb_0, 0))
        self.connect((self.blocks_throttle2_0_0, 0), (self.blocks_add_xx_1_0, 1))
        self.connect((self.blocks_throttle2_0_0, 0), (self.qtgui_sink_x_0_0_0, 0))
        self.connect((self.blocks_unpack_k_bits_bb_1, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.blocks_unpack_k_bits_bb_1, 0), (self.blocks_repeat_0_0_0, 0))
        self.connect((self.blocks_unpack_k_bits_bb_1, 0), (self.blocks_skiphead_0_0_0, 0))
        self.connect((self.digital_chunks_to_symbols_xx_0, 0), (self.blocks_repeat_1, 0))
        self.connect((self.digital_chunks_to_symbols_xx_0, 0), (self.qtgui_time_sink_x_2_0_0_2_0, 0))
        self.connect((self.fec_ber_bf_0_0, 0), (self.qtgui_number_sink_0_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "BFSKloopbacktry")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

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

    def get_m(self):
        return self.m

    def set_m(self, m):
        self.m = m
        self.set_bw((2**self.m)*self.Rb/self.m)
        self.set_fsk_deviation(self.bw/self.m)
        self.set_noise_power(1 / (2 * self.m* 10**(self.eb_n0_dB / 10)))

    def get_Rb(self):
        return self.Rb

    def set_Rb(self, Rb):
        self.Rb = Rb
        self.set_Rs(self.Rb/self.bps)
        self.set_Sps(int((self.bps*self.samp_rate)/self.Rb))
        self.set_bw((2**self.m)*self.Rb/self.m)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_Sps(int((self.bps*self.samp_rate)/self.Rb))
        self.analog_frequency_modulator_fc_0.set_sensitivity(((2*math.pi*self.fsk_deviation)/self.samp_rate))
        self.analog_quadrature_demod_cf_0_0.set_gain((self.samp_rate/(2*math.pi*self.fsk_deviation)))
        self.blocks_multiply_const_vxx_0_0.set_k(self.samp_rate/(2*math.pi*self.fsk_deviation))
        self.blocks_throttle2_0_0.set_sample_rate(self.samp_rate)
        self.qtgui_time_sink_x_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_0_2_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_2_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_2_0_0_2_0.set_samp_rate(self.samp_rate)

    def get_eb_n0_dB(self):
        return self.eb_n0_dB

    def set_eb_n0_dB(self, eb_n0_dB):
        self.eb_n0_dB = eb_n0_dB
        self.set_noise_power(1 / (2 * self.m* 10**(self.eb_n0_dB / 10)))

    def get_bw(self):
        return self.bw

    def set_bw(self, bw):
        self.bw = bw
        self.set_fsk_deviation(self.bw/self.m)
        self.qtgui_sink_x_0_0_0.set_frequency_range(self.center_freq, (self.bw*8))

    def get_bps(self):
        return self.bps

    def set_bps(self, bps):
        self.bps = bps
        self.set_Rs(self.Rb/self.bps)
        self.set_Sps(int((self.bps*self.samp_rate)/self.Rb))

    def get_sig_power(self):
        return self.sig_power

    def set_sig_power(self, sig_power):
        self.sig_power = sig_power

    def get_packet_len(self):
        return self.packet_len

    def set_packet_len(self, packet_len):
        self.packet_len = packet_len
        self.blocks_stream_to_tagged_stream_0.set_packet_len(self.packet_len)
        self.blocks_stream_to_tagged_stream_0.set_packet_len_pmt(self.packet_len)

    def get_num_samples(self):
        return self.num_samples

    def set_num_samples(self, num_samples):
        self.num_samples = num_samples

    def get_noise_power(self):
        return self.noise_power

    def set_noise_power(self, noise_power):
        self.noise_power = noise_power
        self.analog_noise_source_x_0_0.set_amplitude(math.sqrt(self.noise_power))

    def get_ndisp(self):
        return self.ndisp

    def set_ndisp(self, ndisp):
        self.ndisp = ndisp

    def get_fsk_deviation(self):
        return self.fsk_deviation

    def set_fsk_deviation(self, fsk_deviation):
        self.fsk_deviation = fsk_deviation
        self.analog_frequency_modulator_fc_0.set_sensitivity(((2*math.pi*self.fsk_deviation)/self.samp_rate))
        self.analog_quadrature_demod_cf_0_0.set_gain((self.samp_rate/(2*math.pi*self.fsk_deviation)))
        self.blocks_multiply_const_vxx_0_0.set_k(self.samp_rate/(2*math.pi*self.fsk_deviation))

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.qtgui_sink_x_0_0_0.set_frequency_range(self.center_freq, (self.bw*8))

    def get_Sps(self):
        return self.Sps

    def set_Sps(self, Sps):
        self.Sps = Sps
        self.blocks_keep_one_in_n_0.set_n(self.Sps)
        self.blocks_repeat_0_0_0.set_interpolation(self.Sps)
        self.blocks_repeat_0_0_0_0.set_interpolation(self.Sps)
        self.blocks_repeat_1.set_interpolation(self.Sps)

    def get_Rs(self):
        return self.Rs

    def set_Rs(self, Rs):
        self.Rs = Rs




def main(top_block_cls=BFSKloopbacktry, options=None):

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
