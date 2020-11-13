# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 17:18:33 2020

@author: 35192
"""

import re

f = open("seq.txt","r")

def ler_seq (f):
    """
    Parameters
    ----------
    f : TYPE: ficheiro de texto (txt)
        DESCRIPTION: ficheiro de texto aberto 

    Returns
    -------
    seq : TYPE: string
        DESCRIPTION: Função que lê um ficheiro texto aberto linha por linha
                        retornando o seu conteudo
    """
    seq = f.readlines()
    for x in seq:
        print (x)
    return (seq)    

ler_seq(f)        


file = open("sequence.fasta","r")

def ler_FASTA_seq (file):
    """

    Parameters
    ----------
    file : TYPE: FASTA 
        DESCRIPTION: Ficheiro fasta aberto

    Returns
    -------
    sequence : TYPE: string
        DESCRIPTION: Função que lê um ficheiro FASTA linha a linha
                        e retorna o seu conteudo    
    """
    sequence = file.readlines()
    for y in sequence:
        print(y)
    return (sequence)

ler_FASTA_seq(file)
 

#seq = "TTAAATCGCGGTGGAGCATGTG"
seq = input("Digite a sequencia: ")       
def valida (seq):
    """
    Parameters
    ----------
    seq : TYPE: strign
        DESCRIPTION: string que corresponde a uma sequencia de DNA

    Returns
    -------
    bool
        DESCRIPTION: Retorna False se a sequencia conter algum elemento que nao seja os
        nucleotidos A,C,G,T.
        Retorna True se for uma sequencia válida (conter apenas os nucleotidos  A,C,G ou T) 
    """
    seq = seq.upper() 
    valid_dna = "ACGT"
    for letter in seq:
        if letter not in valid_dna:  
            return False
    return True
              
valida (seq)    



#seq = "TTAAATCGCGGTGGAGCATGTG"
#seq = "ACFG"
seq = input("Digite a sequencia: ")  

def complemento_inverso (seq):
    """
    Parameters
    ----------
    seq : TYPE: string
        DESCRIPTION: string com a sequencia de DNA que 
                      se pretende obter o complemento inverso

    Returns
    -------
    seq_comp_inv : TYPE: string
        DESCRIPTION: string que retorna o complemento inverso da sequencia inserida
    """ 
    if valida(seq) is False:
            raise TypeError("A sequencia só pode conter ACGT")    
    else:
        seq_comp_inv = seq[::-1].lower()
        seq_comp_inv = seq_comp_inv.replace('a','T')
        seq_comp_inv = seq_comp_inv.replace('t','A')
        seq_comp_inv = seq_comp_inv.replace('g','C')
        seq_comp_inv = seq_comp_inv.replace('c','G')
    return (seq_comp_inv)

complemento_inverso(seq)


seq = input("Digite a sequencia: ")  

def transcricao(seq):
    """
    Parameters
    ----------
    seq : TYPE: string 
        DESCRIPTION:string com a sequencia de DNA que 
                      se pretende obter a transcricao

    Returns
    -------
    seq_trans : TYPE: string
        DESCRIPTION: transcrito da sequencia de DNA
        
    Funçao que recebe uma string (seq) e retorna outra string(seq_trans)   
    Esta funçao realiza a transcriçao, que é a substituiçao dos nucleotidos T por U 
    """
    if valida(seq) is False:
        raise TypeError("A sequencia só pode conter ACGT")  
    else:    
        seq_trans = seq.upper().replace('T','U')
    return(seq_trans)

transcricao(seq)


#seq = "ACFG"
seq = input("Digite a sequencia: ")  

def traducao (seq):
    """
    Parameters
    ----------
    seq : TYPE: string
        DESCRIPTION: string com a sequencia de DNA que se
                     pretende traduzir em aminoácidos

    Raises
    ------
    TypeError
        DESCRIPTION: Se nao for uma sequencia válida de DNA (nao conter apenas A,C,G,T)

    Returns
    -------
    cadeia_aminoacidos : TYPE: string
        DESCRIPTION: sequencia de aminoacidos
        
   Funçao que atraves de uma biblioteca (gencode) atribui a cada codao o aminoacido
   correspondente. Retornando uma sequencia de aminoacidos (cadeia_aminoacidos)     
   
    """ 
    if (len(seq) < 3 ):
        raise TypeError("A sequencia nao pode ser inferior a 3")
    else:    
        if valida(seq) is False:
            raise TypeError("A sequencia só pode conter ACGT")     
        elif valida(seq) is True:    
            cadeia_DNA_separada = [(seq[i:i + 3]) for i in range(0, len(seq), 3)]
            gencode = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
            
            cadeia_aminoacidos = ""
            for cadeia in cadeia_DNA_separada:
                cadeia_aminoacidos += gencode.get(cadeia,"")
            return cadeia_aminoacidos

traducao(seq)
        

seq = input("Digite a sequencia ")

def contar_bases (seq):
    """
    Parameters
    ----------
    seq : TYPE: string
        DESCRIPTION: sequencia de DNA,RNA ou aminoacidos 

    Returns
    -------
    bases : TYPE: dicionario 
        DESCRIPTION: retorna um dicionario com a contagem de cada elemento da sequencia
   Funçao que recebe uma sequencia de DNA,RNA ou aminoacidos 
   e retorna a contagem de cada elemento num dicionario.     
    """
    A, C, G, T, U = 0, 0, 0, 0, 0
    seq = seq.upper()
    if "ACGTU" not in seq:   
        raise TypeError("A sequencia nao e valida") #NAO FUNCIONA
    else:    
        if "T" in seq:
            A = seq.count('A')
            C = seq.count('C')
            G = seq.count('G')
            T = seq.count('T')
            bases = {"A": A , "C": C , "G": G, "T": T}
        elif "U" in seq:
            A = seq.count('A')
            C = seq.count('C')
            G = seq.count('G')
            U = seq.count('U')
            bases = {"A": A , "C": C , "G": G, "U": U}        
        elif "A" or "C" or "G" in seq:
            A = seq.count('A')
            C = seq.count('C')
            G = seq.count('G')
            bases = {"A": A , "C": C , "G": G}
    
    return bases

contar_bases(seq)    

#sequencia = "CGATGGCTGATTGAATGGCGTAAATGTAACATGGCTTAAAAATGGCGTGGTAAAATGGACGCATGACTATGTAAATGTAAATGGCGTACTAGTATGGCTAGGTACTAGCGATGTAG"
sequencia = input("Digite a sequencia ")

def traducao_RNA (sequencia):
    """
    Parameters
    ----------
    sequencia : TYPE: string
        DESCRIPTION: string que recebe uma sequencia de RNA

    Raises
    ------
    TypeError
        DESCRIPTION: Se nao for uma sequencia válida de RNA (nao conter apenas A,C,G ou U)

    Returns
    -------
    cadeia_aminoacidos : TYPE string
        DESCRIPTION:
    Funçao que atraves de uma biblioteca (gencodeRNA) atribui a cada codao o aminoacido
    correspondente. Retornando uma sequencia de aminoacidos (cadeia_aminoacidos)           

    """
    sequencia = sequencia.upper() 
    valid_rna = "ACGU"
    if (len(sequencia) < 3 ):
        raise TypeError("A sequencia nao pode ser inferior a 3")
    else:    
        for letter in sequencia:
            if letter not in valid_rna:
                raise TypeError("A sequencia só pode conter A,C,G ou U")            
    cadeia_RNA_separada = [(sequencia[i:i + 3]) for i in range(0, len(sequencia), 3)]
    gencodeRNA = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W'}
    cadeia_aminoacidos = ""
    for cadeia in cadeia_RNA_separada:
        cadeia_aminoacidos += gencodeRNA.get(cadeia,"")
    return cadeia_aminoacidos
    
traducao_RNA(sequencia)    
    
    
sequencia = input("Digite a sequencia ")

def reading_frames (sequencia):
    """
    Parameters
    ----------
    sequencia : TYPE: string
        DESCRIPTION: string que contem uma sequencia de DNA 

    Returns
    -------
    reading_frames : TYPE: lista
        DESCRIPTION: Retorna uma Lista com as reading frames
    Função que recebe uma sequencia de DNA e 
    devolve uma lista com as reding frames    

    """
    if valida(sequencia) is True: # So esta a validar uma sequencia de DNA
        reading_frames = re.findall('...', sequencia) 
    else:
        raise TypeError("A sequencia é inválida")
    return (reading_frames)
   
reading_frames(sequencia)    


sequencia = input("Digite a sequencia: ")

def get_proteins(sequencia):
    """
    Parameters
    ----------
    sequencia : TYPE: string
        DESCRIPTION. string de uma sequencia de DNA
        
    Returns
    -------
    palavras: TYPE: lista
        DESCRIPTION: Retorna uma lista de proteinas 
   Funçao que recebe uma sequencia de DNA
   e devolve a lista de todas as proteinas ordenadas por tamanho
   e ordem alfabetica para as do mesmo tamanho
    """
    palavras_encontradas = []
    frames = reading_frames(sequencia)
    lista = "".join(frames) # converter a lista para uma strign pois a funçao traducao recebe uma string
    tradu = []
    tradu.append(traducao(lista))
    for palavra in tradu: # tradu e uma lista onde estao os aminoacidos
        caracteres = [char for char in palavra]
        palavra_encontrada = ""                 
        escreve = False                         
        for caracter in caracteres:    # percorre a lista tradu         
            if caracter == "M":
                escreve = True
            if caracter == "_" and len(palavra_encontrada) != 0: 
                escreve = False                                 
                palavra_encontrada += caracter                  
                palavras_encontradas.append(palavra_encontrada) 
                palavra_encontrada = ""                              
            if escreve:
                palavra_encontrada += caracter                    
    palavras_encontradas_ord = sorted(palavras_encontradas) 
    palavras_s_duplicados = sorted(set(palavras_encontradas_ord), key= len, reverse = True)
    
    return(palavras_s_duplicados)   
    
get_proteins(sequencia)    
    
    

















