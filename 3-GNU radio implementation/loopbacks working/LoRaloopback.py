#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: loopback LoRa
# Author: clsor
# GNU Radio version: 3.10.10.0
"""This script provides the blocks for the LoRa loopback modulation-demodulation with adding of AWGN."""
from PyQt5 import Qt
from gnuradio import qtgui
from PyQt5 import QtCore
from gnuradio import analog
from gnuradio import blocks
import numpy
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
import gnuradio.lora_sdr as lora_sdr
import math
import sip
import untitled_epy_block_0 as epy_block_0  # embedded python block



class untitled(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "loopback LoRa", catch_exceptions=True)
        Qt.QWidget.__init__(self)
        self.setWindowTitle("loopback LoRa")
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

        self.settings = Qt.QSettings("GNU Radio", "untitled")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

   # Variables
        self.sf = sf = 7
        self.bw = bw = 125000
        self.sig_power = sig_power = 1
        self.packet_len = packet_len = 16*15
        self.eb_n0_dB = eb_n0_dB = 15
        self.Rb = Rb = (sf*bw)/(2**sf)
        self.soft_decoding = soft_decoding = False
        self.samp_rate = samp_rate = bw
        self.preamb_len = preamb_len = 8
        self.pay_len = pay_len = packet_len
        self.num_samples = num_samples = 100000
        self.noise_power = noise_power = ((sig_power*bw)/Rb)*10**(-eb_n0_dB/10)
        self.ndisp = ndisp = 32
        self.impl_head = impl_head = False
        self.has_crc = has_crc = False
        self.cr = cr = 0
        self.clk_offset = clk_offset = 0
        self.center_freq = center_freq = 868100000
        self.Sps = Sps = 40
        self.Rs = Rs = Rb/sf
        self.Idro = Idro = False

       # Blocks
        self.qtgui_time_sink_x_2_0 = qtgui.time_sink_f(
            ndisp, #size
            samp_rate, #samp_rate
            'Comparison bits', #name
            2, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_2_0.set_update_time(0.0010)
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


        for i in range(2):
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
        self.qtgui_time_sink_x_0_2_2_0_0 = qtgui.time_sink_c(
            ndisp, #size
            samp_rate, #samp_rate
            'modulate', #name
            1, #number of inputs
            None # parent
        )
        self.qtgui_time_sink_x_0_2_2_0_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0_2_2_0_0.set_y_axis(-0.5, 1.5)

        self.qtgui_time_sink_x_0_2_2_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0_2_2_0_0.enable_tags(True)
        self.qtgui_time_sink_x_0_2_2_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0_2_2_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0_2_2_0_0.enable_grid(True)
        self.qtgui_time_sink_x_0_2_2_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0_2_2_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0_2_2_0_0.enable_stem_plot(False)

        self.qtgui_time_sink_x_0_2_2_0_0.disable_legend()

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


        for i in range(2):
            if len(labels[i]) == 0:
                if (i % 2 == 0):
                    self.qtgui_time_sink_x_0_2_2_0_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_0_2_2_0_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_0_2_2_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_2_2_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_2_2_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_2_2_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_2_2_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_2_2_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_2_2_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0_2_2_0_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_2_2_0_0_win)
        self.qtgui_sink_x_0_0 = qtgui.sink_c(
            1024, #fftsize
            window.WIN_BLACKMAN_hARRIS, #wintype
            center_freq, #fc
            (samp_rate*8), #bw
            'before noise', #name
            True, #plotfreq
            True, #plotwaterfall
            True, #plottime
            True, #plotconst
            None # parent
        )
        self.qtgui_sink_x_0_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_0_win = sip.wrapinstance(self.qtgui_sink_x_0_0.qwidget(), Qt.QWidget)

        self.qtgui_sink_x_0_0.enable_rf_freq(True)

        self.top_layout.addWidget(self._qtgui_sink_x_0_0_win)
        self.qtgui_sink_x_0 = qtgui.sink_c(
            1024, #fftsize
            window.WIN_BLACKMAN_hARRIS, #wintype
            center_freq, #fc
            (samp_rate*8), #bw
            'after noise', #name
            True, #plotfreq
            True, #plotwaterfall
            True, #plottime
            True, #plotconst
            None # parent
        )
        self.qtgui_sink_x_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_win = sip.wrapinstance(self.qtgui_sink_x_0.qwidget(), Qt.QWidget)

        self.qtgui_sink_x_0.enable_rf_freq(False)

        self.top_layout.addWidget(self._qtgui_sink_x_0_win)
        self.qtgui_number_sink_0 = qtgui.number_sink(
            gr.sizeof_float,
            0,
            qtgui.NUM_GRAPH_HORIZ,
            1,
            None # parent
        )
        self.qtgui_number_sink_0.set_update_time(0.10)
        self.qtgui_number_sink_0.set_title('BER between Tx and Rx')

        labels = ['', '', '', '', '',
            '', '', '', '', '']
        units = ['', '', '', '', '',
            '', '', '', '', '']
        colors = [("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"),
            ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black")]
        factor = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]

        for i in range(1):
            self.qtgui_number_sink_0.set_min(i, -7)
            self.qtgui_number_sink_0.set_max(i, 1)
            self.qtgui_number_sink_0.set_color(i, colors[i][0], colors[i][1])
            if len(labels[i]) == 0:
                self.qtgui_number_sink_0.set_label(i, "Data {0}".format(i))
            else:
                self.qtgui_number_sink_0.set_label(i, labels[i])
            self.qtgui_number_sink_0.set_unit(i, units[i])
            self.qtgui_number_sink_0.set_factor(i, factor[i])

        self.qtgui_number_sink_0.enable_autoscale(False)
        self._qtgui_number_sink_0_win = sip.wrapinstance(self.qtgui_number_sink_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_number_sink_0_win)
        self.lora_sdr_whitening_0 = lora_sdr.whitening(False,True,';','packet_len')
        self.lora_sdr_modulate_0 = lora_sdr.modulate(sf, samp_rate, bw, [0x12], (int(20*2**sf*samp_rate/bw)),preamb_len)
        self.lora_sdr_interleaver_0 = lora_sdr.interleaver(cr, sf, 0, 125000)
        self.lora_sdr_header_decoder_0 = lora_sdr.header_decoder(impl_head, 3, 255, False, 0, True)
        self.lora_sdr_header_0 = lora_sdr.header(impl_head, has_crc, cr)
        self.lora_sdr_hamming_enc_0 = lora_sdr.hamming_enc(cr, sf)
        self.lora_sdr_hamming_dec_0 = lora_sdr.hamming_dec(soft_decoding)
        self.lora_sdr_gray_mapping_0 = lora_sdr.gray_mapping( soft_decoding)
        self.lora_sdr_gray_demap_0 = lora_sdr.gray_demap(sf)
        self.lora_sdr_frame_sync_0 = lora_sdr.frame_sync(int(center_freq), bw, sf, impl_head, [18], (int(samp_rate/bw)),preamb_len)
        self.lora_sdr_fft_demod_1 = lora_sdr.fft_demod( soft_decoding, True)
        self.lora_sdr_dewhitening_0 = lora_sdr.dewhitening()
        self.lora_sdr_deinterleaver_0 = lora_sdr.deinterleaver( soft_decoding)
        self.lora_sdr_add_crc_0 = lora_sdr.add_crc(has_crc)
        self.fec_ber_bf_0 = fec.ber_bf(False, 100, -7.0)
        self._eb_n0_dB_range = qtgui.Range(-5, 15, 1/100, 15, 200)
        self._eb_n0_dB_win = qtgui.RangeWidget(self._eb_n0_dB_range, self.set_eb_n0_dB, "Eb/N0", "counter_slider", float, QtCore.Qt.Horizontal)
        self.top_grid_layout.addWidget(self._eb_n0_dB_win, 2, 0, 1, 1)
        for r in range(2, 3):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.blocks_throttle2_1 = blocks.throttle( gr.sizeof_gr_complex*1, samp_rate, True, 0 if "auto" == "auto" else max( int(float(0.1) * samp_rate) if "auto" == "time" else int(0.1), 1) )
        self.blocks_threshold_ff_0 = blocks.threshold_ff(0.5, 0.5, 0)
        self.blocks_stream_to_tagged_stream_2 = blocks.stream_to_tagged_stream(gr.sizeof_char, 1, packet_len, "packet_len")
        self.blocks_skiphead_0_0 = blocks.skiphead(gr.sizeof_char*1, 0)
        self.blocks_skiphead_0 = blocks.skiphead(gr.sizeof_char*1, 0)
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_cc(1)
        self.blocks_float_to_char_0 = blocks.float_to_char(1, 1)
        self.blocks_delay_0 = blocks.delay(gr.sizeof_char*1, 0)
        self.blocks_char_to_float_0_1 = blocks.char_to_float(1, 1)
        self.blocks_char_to_float_0_0_0_0 = blocks.char_to_float(1, 1)
        self.blocks_char_to_float_0_0_0 = blocks.char_to_float(1, 1)
        self.blocks_add_xx_1 = blocks.add_vcc(1)
        self.analog_random_source_x_0 = blocks.vector_source_b(list(map(int, numpy.random.randint(0, 2, num_samples))), False)
        self.analog_noise_source_x_0_0 = analog.noise_source_c(analog.GR_GAUSSIAN, (math.sqrt(2*noise_power)), 0)


        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.lora_sdr_header_decoder_0, 'frame_info'), (self.lora_sdr_frame_sync_0, 'frame_info'))
        self.connect((self.analog_noise_source_x_0_0, 0), (self.blocks_multiply_const_vxx_0, 0))
        self.connect((self.analog_random_source_x_0, 0), (self.blocks_skiphead_0_0, 0))
        self.connect((self.analog_random_source_x_0, 0), (self.blocks_stream_to_tagged_stream_2, 0))
        self.connect((self.blocks_add_xx_1, 0), (self.lora_sdr_frame_sync_0, 0))
        self.connect((self.blocks_add_xx_1, 0), (self.qtgui_sink_x_0, 0))
        self.connect((self.blocks_char_to_float_0_0_0, 0), (self.qtgui_time_sink_x_2_0, 1))
        self.connect((self.blocks_char_to_float_0_0_0_0, 0), (self.blocks_threshold_ff_0, 0))
        self.connect((self.blocks_char_to_float_0_1, 0), (self.qtgui_time_sink_x_2_0, 0))
        self.connect((self.blocks_delay_0, 0), (self.blocks_char_to_float_0_1, 0))
        self.connect((self.blocks_delay_0, 0), (self.fec_ber_bf_0, 0))
        self.connect((self.blocks_float_to_char_0, 0), (self.blocks_skiphead_0, 0))
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.blocks_add_xx_1, 0))
        self.connect((self.blocks_skiphead_0, 0), (self.blocks_char_to_float_0_0_0, 0))
        self.connect((self.blocks_skiphead_0, 0), (self.fec_ber_bf_0, 1))
        self.connect((self.blocks_skiphead_0_0, 0), (self.blocks_delay_0, 0))
        self.connect((self.blocks_stream_to_tagged_stream_2, 0), (self.lora_sdr_whitening_0, 0))
        self.connect((self.blocks_threshold_ff_0, 0), (self.blocks_float_to_char_0, 0))
        self.connect((self.blocks_throttle2_1, 0), (self.blocks_add_xx_1, 1))
        self.connect((self.blocks_throttle2_1, 0), (self.qtgui_sink_x_0_0, 0))
        self.connect((self.fec_ber_bf_0, 0), (self.qtgui_number_sink_0, 0))
        self.connect((self.lora_sdr_add_crc_0, 0), (self.lora_sdr_hamming_enc_0, 0))
        self.connect((self.lora_sdr_deinterleaver_0, 0), (self.lora_sdr_hamming_dec_0, 0))
        self.connect((self.lora_sdr_dewhitening_0, 0), (self.blocks_char_to_float_0_0_0_0, 0))
        self.connect((self.lora_sdr_fft_demod_1, 0), (self.lora_sdr_gray_mapping_0, 0))
        self.connect((self.lora_sdr_frame_sync_0, 0), (self.lora_sdr_fft_demod_1, 0))
        self.connect((self.lora_sdr_gray_demap_0, 0), (self.lora_sdr_modulate_0, 0))
        self.connect((self.lora_sdr_gray_mapping_0, 0), (self.lora_sdr_deinterleaver_0, 0))
        self.connect((self.lora_sdr_hamming_dec_0, 0), (self.lora_sdr_header_decoder_0, 0))
        self.connect((self.lora_sdr_hamming_enc_0, 0), (self.lora_sdr_interleaver_0, 0))
        self.connect((self.lora_sdr_header_0, 0), (self.lora_sdr_add_crc_0, 0))
        self.connect((self.lora_sdr_header_decoder_0, 0), (self.lora_sdr_dewhitening_0, 0))
        self.connect((self.lora_sdr_interleaver_0, 0), (self.lora_sdr_gray_demap_0, 0))
        self.connect((self.lora_sdr_modulate_0, 0), (self.blocks_throttle2_1, 0))
        self.connect((self.lora_sdr_modulate_0, 0), (self.qtgui_time_sink_x_0_2_2_0_0, 0))
        self.connect((self.lora_sdr_whitening_0, 0), (self.lora_sdr_header_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "untitled")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_sf(self):
        return self.sf

    def set_sf(self, sf):
        self.sf = sf
        self.set_Rb((self.sf*self.bw)/(2**self.sf))
        self.set_Rs(self.Rb/self.sf)
        self.lora_sdr_gray_demap_0.set_sf(self.sf)
        self.lora_sdr_hamming_enc_0.set_sf(self.sf)
        self.lora_sdr_interleaver_0.set_sf(self.sf)
        self.lora_sdr_modulate_0.set_sf(self.sf)

    def get_bw(self):
        return self.bw

    def set_bw(self, bw):
        self.bw = bw
        self.set_Rb((self.sf*self.bw)/(2**self.sf))
        self.set_noise_power(((self.sig_power*self.bw)/self.Rb)*10**(-self.eb_n0_dB/10))
        self.set_samp_rate(self.bw)

    def get_sig_power(self):
        return self.sig_power

    def set_sig_power(self, sig_power):
        self.sig_power = sig_power
        self.set_noise_power(((self.sig_power*self.bw)/self.Rb)*10**(-self.eb_n0_dB/10))

    def get_packet_len(self):
        return self.packet_len

    def set_packet_len(self, packet_len):
        self.packet_len = packet_len
        self.set_pay_len(self.packet_len)
        self.blocks_stream_to_tagged_stream_2.set_packet_len(self.packet_len)
        self.blocks_stream_to_tagged_stream_2.set_packet_len_pmt(self.packet_len)

    def get_eb_n0_dB(self):
        return self.eb_n0_dB

    def set_eb_n0_dB(self, eb_n0_dB):
        self.eb_n0_dB = eb_n0_dB
        self.set_noise_power(((self.sig_power*self.bw)/self.Rb)*10**(-self.eb_n0_dB/10))

    def get_Rb(self):
        return self.Rb

    def set_Rb(self, Rb):
        self.Rb = Rb
        self.set_Rs(self.Rb/self.sf)
        self.set_noise_power(((self.sig_power*self.bw)/self.Rb)*10**(-self.eb_n0_dB/10))

    def get_soft_decoding(self):
        return self.soft_decoding

    def set_soft_decoding(self, soft_decoding):
        self.soft_decoding = soft_decoding

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.blocks_throttle2_1.set_sample_rate(self.samp_rate)
        self.qtgui_sink_x_0.set_frequency_range(self.center_freq, (self.samp_rate*8))
        self.qtgui_sink_x_0_0.set_frequency_range(self.center_freq, (self.samp_rate*8))
        self.qtgui_time_sink_x_0_2_2_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_2_0.set_samp_rate(self.samp_rate)

    def get_preamb_len(self):
        return self.preamb_len

    def set_preamb_len(self, preamb_len):
        self.preamb_len = preamb_len

    def get_pay_len(self):
        return self.pay_len

    def set_pay_len(self, pay_len):
        self.pay_len = pay_len

    def get_num_samples(self):
        return self.num_samples

    def set_num_samples(self, num_samples):
        self.num_samples = num_samples

    def get_noise_power(self):
        return self.noise_power

    def set_noise_power(self, noise_power):
        self.noise_power = noise_power
        self.analog_noise_source_x_0_0.set_amplitude((math.sqrt(2*self.noise_power)))

    def get_ndisp(self):
        return self.ndisp

    def set_ndisp(self, ndisp):
        self.ndisp = ndisp

    def get_impl_head(self):
        return self.impl_head

    def set_impl_head(self, impl_head):
        self.impl_head = impl_head

    def get_has_crc(self):
        return self.has_crc

    def set_has_crc(self, has_crc):
        self.has_crc = has_crc

    def get_cr(self):
        return self.cr

    def set_cr(self, cr):
        self.cr = cr
        self.lora_sdr_hamming_enc_0.set_cr(self.cr)
        self.lora_sdr_header_0.set_cr(self.cr)
        self.lora_sdr_interleaver_0.set_cr(self.cr)

    def get_clk_offset(self):
        return self.clk_offset

    def set_clk_offset(self, clk_offset):
        self.clk_offset = clk_offset

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.qtgui_sink_x_0.set_frequency_range(self.center_freq, (self.samp_rate*8))
        self.qtgui_sink_x_0_0.set_frequency_range(self.center_freq, (self.samp_rate*8))

    def get_Sps(self):
        return self.Sps

    def set_Sps(self, Sps):
        self.Sps = Sps

    def get_Rs(self):
        return self.Rs

    def set_Rs(self, Rs):
        self.Rs = Rs

    def get_Idro(self):
        return self.Idro

    def set_Idro(self, Idro):
        self.Idro = Idro




def main(top_block_cls=untitled, options=None):

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
