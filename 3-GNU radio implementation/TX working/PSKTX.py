#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: PSKTX
# Author: clsor
# GNU Radio version: 3.10.10.0

from PyQt5 import Qt
from gnuradio import qtgui
from PyQt5 import QtCore
from gnuradio import analog
from gnuradio import blocks
import pmt
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
import PSKTX_wform as wform  # embedded python module
import math
import sip



class PSKTX(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "PSKTX", catch_exceptions=True)
        Qt.QWidget.__init__(self)
        self.setWindowTitle("PSKTX")
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

        self.settings = Qt.QSettings("GNU Radio", "PSKTX")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Variables
        ##################################################
        self.constellation = constellation = digital.constellation_qpsk().points()
        self.sf_lora = sf_lora = 7
        self.bw_lora = bw_lora = 125000
        self.M = M = len(constellation )
        self.bps = bps = int(math.log(M,2))
        self.Rb = Rb = (sf_lora*bw_lora)/2**sf_lora
        self.m = m = int(math.log2(M))
        self.Sps = Sps = 40
        self.Rs = Rs = Rb/bps
        self.sig_power = sig_power = 1
        self.samp_rate = samp_rate = Rs*Sps
        self.ntaps = ntaps = 16*Sps
        self.eb_n0_dB = eb_n0_dB = 15
        self.bw = bw = (2**m)*Rb/m
        self.beta = beta = 1
        self.packet_len = packet_len = 240
        self.num_samples = num_samples = 100000
        self.noise_power = noise_power = (sig_power*bw/(Sps*Rs))*10**(-eb_n0_dB/10)
        self.h = h = wform.rrcos(Sps,ntaps,beta)
        self.delay = delay = 16
        self.center_freq = center_freq = 100000
        self.Fmax = Fmax = samp_rate/2

        ##################################################
        # Blocks
        ##################################################

        self.soapy_hackrf_sink_0 = None
        dev = 'driver=hackrf'
        stream_args = ''
        tune_args = ['']
        settings = ['']

        self.soapy_hackrf_sink_0 = soapy.sink(dev, "fc32", 1, 'driver=hackrf',
                                  stream_args, tune_args, settings)
        self.soapy_hackrf_sink_0.set_sample_rate(0, samp_rate)
        self.soapy_hackrf_sink_0.set_bandwidth(0, bw)
        self.soapy_hackrf_sink_0.set_frequency(0, center_freq)
        self.soapy_hackrf_sink_0.set_gain(0, 'AMP', False)
        self.soapy_hackrf_sink_0.set_gain(0, 'VGA', min(max(16, 0.0), 47.0))
        self.qtgui_sink_x_0_0 = qtgui.sink_c(
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
        self.interp_fir_filter_xxx_0 = filter.interp_fir_filter_ccf(Sps, h)
        self.interp_fir_filter_xxx_0.declare_sample_delay(0)
        self._eb_n0_dB_range = qtgui.Range(-5, 15, 1/100, 15, 200)
        self._eb_n0_dB_win = qtgui.RangeWidget(self._eb_n0_dB_range, self.set_eb_n0_dB, "Eb/N0", "counter_slider", float, QtCore.Qt.Horizontal)
        self.top_grid_layout.addWidget(self._eb_n0_dB_win, 2, 0, 1, 1)
        for r in range(2, 3):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.digital_chunks_to_symbols_xx_0 = digital.chunks_to_symbols_bc(constellation, 1)
        self.blocks_throttle2_0 = blocks.throttle( gr.sizeof_gr_complex*1, samp_rate, True, 0 if "auto" == "auto" else max( int(float(0.1) * samp_rate) if "auto" == "time" else int(0.1), 1) )
        self.blocks_stream_to_tagged_stream_0 = blocks.stream_to_tagged_stream(gr.sizeof_char, 1, packet_len, "packet_len")
        self.blocks_packed_to_unpacked_xx_0 = blocks.packed_to_unpacked_bb(bps, gr.GR_MSB_FIRST)
        self.blocks_pack_k_bits_bb_0 = blocks.pack_k_bits_bb(8)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_char*1, 'C:\\Users\\clsor\\OneDrive\\Documents\\MATLAB\\Master_Thesis_Clara\\Master_Thesis_Clara\\3-GNU radio implementation\\SDR files of bits\\sdrinput.bin', False, 0, 0)
        self.blocks_file_source_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_add_xx_0 = blocks.add_vcc(1)
        self.analog_noise_source_x_0 = analog.noise_source_c(analog.GR_GAUSSIAN, (math.sqrt(2*noise_power)), 0)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_noise_source_x_0, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.qtgui_sink_x_0, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.soapy_hackrf_sink_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.blocks_stream_to_tagged_stream_0, 0))
        self.connect((self.blocks_pack_k_bits_bb_0, 0), (self.blocks_packed_to_unpacked_xx_0, 0))
        self.connect((self.blocks_packed_to_unpacked_xx_0, 0), (self.digital_chunks_to_symbols_xx_0, 0))
        self.connect((self.blocks_stream_to_tagged_stream_0, 0), (self.blocks_pack_k_bits_bb_0, 0))
        self.connect((self.blocks_throttle2_0, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.blocks_throttle2_0, 0), (self.qtgui_sink_x_0_0, 0))
        self.connect((self.digital_chunks_to_symbols_xx_0, 0), (self.interp_fir_filter_xxx_0, 0))
        self.connect((self.interp_fir_filter_xxx_0, 0), (self.blocks_throttle2_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "PSKTX")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_constellation(self):
        return self.constellation

    def set_constellation(self, constellation):
        self.constellation = constellation
        self.set_M(len(self.constellation ))
        self.digital_chunks_to_symbols_xx_0.set_symbol_table(self.constellation)

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

    def get_bps(self):
        return self.bps

    def set_bps(self, bps):
        self.bps = bps
        self.set_Rs(self.Rb/self.bps)

    def get_Rb(self):
        return self.Rb

    def set_Rb(self, Rb):
        self.Rb = Rb
        self.set_Rs(self.Rb/self.bps)
        self.set_bw((2**self.m)*self.Rb/self.m)

    def get_m(self):
        return self.m

    def set_m(self, m):
        self.m = m
        self.set_bw((2**self.m)*self.Rb/self.m)

    def get_Sps(self):
        return self.Sps

    def set_Sps(self, Sps):
        self.Sps = Sps
        self.set_h(wform.rrcos(self.Sps,self.ntaps,self.beta))
        self.set_noise_power((self.sig_power*self.bw/(self.Sps*self.Rs))*10**(-self.eb_n0_dB/10))
        self.set_ntaps(16*self.Sps)
        self.set_samp_rate(self.Rs*self.Sps)

    def get_Rs(self):
        return self.Rs

    def set_Rs(self, Rs):
        self.Rs = Rs
        self.set_noise_power((self.sig_power*self.bw/(self.Sps*self.Rs))*10**(-self.eb_n0_dB/10))
        self.set_samp_rate(self.Rs*self.Sps)

    def get_sig_power(self):
        return self.sig_power

    def set_sig_power(self, sig_power):
        self.sig_power = sig_power
        self.set_noise_power((self.sig_power*self.bw/(self.Sps*self.Rs))*10**(-self.eb_n0_dB/10))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_Fmax(self.samp_rate/2)
        self.blocks_throttle2_0.set_sample_rate(self.samp_rate)
        self.qtgui_sink_x_0.set_frequency_range(self.center_freq, (self.samp_rate*8))
        self.soapy_hackrf_sink_0.set_sample_rate(0, self.samp_rate)

    def get_ntaps(self):
        return self.ntaps

    def set_ntaps(self, ntaps):
        self.ntaps = ntaps
        self.set_h(wform.rrcos(self.Sps,self.ntaps,self.beta))

    def get_eb_n0_dB(self):
        return self.eb_n0_dB

    def set_eb_n0_dB(self, eb_n0_dB):
        self.eb_n0_dB = eb_n0_dB
        self.set_noise_power((self.sig_power*self.bw/(self.Sps*self.Rs))*10**(-self.eb_n0_dB/10))

    def get_bw(self):
        return self.bw

    def set_bw(self, bw):
        self.bw = bw
        self.set_noise_power((self.sig_power*self.bw/(self.Sps*self.Rs))*10**(-self.eb_n0_dB/10))
        self.qtgui_sink_x_0_0.set_frequency_range(self.center_freq, (self.bw*8))
        self.soapy_hackrf_sink_0.set_bandwidth(0, self.bw)

    def get_beta(self):
        return self.beta

    def set_beta(self, beta):
        self.beta = beta
        self.set_h(wform.rrcos(self.Sps,self.ntaps,self.beta))

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
        self.analog_noise_source_x_0.set_amplitude((math.sqrt(2*self.noise_power)))

    def get_h(self):
        return self.h

    def set_h(self, h):
        self.h = h
        self.interp_fir_filter_xxx_0.set_taps(self.h)

    def get_delay(self):
        return self.delay

    def set_delay(self, delay):
        self.delay = delay

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.qtgui_sink_x_0.set_frequency_range(self.center_freq, (self.samp_rate*8))
        self.qtgui_sink_x_0_0.set_frequency_range(self.center_freq, (self.bw*8))
        self.soapy_hackrf_sink_0.set_frequency(0, self.center_freq)

    def get_Fmax(self):
        return self.Fmax

    def set_Fmax(self, Fmax):
        self.Fmax = Fmax




def main(top_block_cls=PSKTX, options=None):

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
