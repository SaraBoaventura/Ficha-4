# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 10:54:18 2020

@author: 35192
"""


import unittest

import Ficha4


class TestFicha4 (unittest.TestCase):
    def test_ler_seq (self):
       self.assertEqual(Ficha4.ler_seq("ACGT"),'ACGT')
    def test_ler_fasta_seq(self):
        self.assertEqual(Ficha4.ler_FASTA_seq("ACGT"),'ACGT')
    def test_complemento_inverso(self):
        self.assertEqual(Ficha4.complemento_inverso("AATCGATCG"),'CGATCGATT')
        self.assertRaises(TypeError, Ficha4.complemento_inverso, 'AACGTE')
        self.assertRaises(TypeError, Ficha4.complemento_inverso, 'AACGT1')
        self.assertRaises(TypeError, Ficha4.complemento_inverso, 'AACGT?')
    def test_transcricao(self):
        self.assertEqual(Ficha4.transcricao("AATCGATCG"),'AAUCGAUCG')
        self.assertRaises(TypeError, Ficha4.transcricao, 'AACGTE')
        self.assertRaises(TypeError, Ficha4.transcricao, 'AACGT1')
        self.assertRaises(TypeError, Ficha4.transcricao, 'AACGT?')
    def test_traducao(self):
        self.assertEqual(Ficha4.traducao(""),'')
        self.assertRaises(TypeError, Ficha4.traducao, 'AACGTE')
        self.assertRaises(TypeError, Ficha4.traducao, 'AACGT1')
        self.assertRaises(TypeError, Ficha4.traducao, 'AACGT?')
    def test_valida(self):
        self.assertEqual(Ficha4.valida("ACGTCCGT"), True)
        self.assertEqual(Ficha4.valida("ACGTECGT"), False)
        self.assertEqual(Ficha4.valida("ACGTCGT1"), False)
    def test_contar_bases(self):
        self.assertRaises(TypeError, Ficha4.contar_bases, 'ACFD')
    def test_traducao_RNA(self):
        self.assertEqual(Ficha4.traducao_RNA("AUC"), 'I')
        self.assertRaises(TypeError, Ficha4.traducao_RNA, 'AU')
        self.assertRaises(TypeError, Ficha4.traducao_RNA, 'ACGE')
        self.assertRaises(TypeError, Ficha4.traducao_RNA, 'ACGT')
    def test_reading_frames(self):
        self.assertRaises(TypeError, Ficha4.reading_frames, 'ACGE')
    def test_get_proteins(self):
        self.assertEqual(Ficha4.get_proteins("AUC"), 'I')   # ???????

                   
if __name__ == '__main__':
    unittest.main()
    
    
    
    
    
    
    
    
        