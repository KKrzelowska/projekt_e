from przydatne_funkcje import Score, profil, podobienstwo
import random as rand

def randomized_motif_search(DNA, k, t):
    """Funkcja zwracająca motyw otrzymany poprzez zastosowanie alorytmu losowego
       k - dlugosc meru
       DNA - tabelka sekwencji nukleotydowych
       t - ilosc sekwencji w dna """

    Motifs = []
    for x in range(0, t):                                      #Petla iterujaca po ilosci sekwencji w dna
        Motifs.append(rand.randint(0, len(DNA[0])) - k)        #Nadawanie losowych wartosci indeksow
    bestMotifs = Motifs                                        #Ktore nadajemy rownież zmiennej bestMotifs
    while 1:                                                   #Nieskonczona pętla
        macierz = profil(Motifs, DNA, k)                       #Zwracanie profilu ze zmiennej Motifs
        for x in range(0, len(Motifs)):                        #Petla od zera do ilosci elementow zmiennej Motifs
            Motifs[x] = podobienstwo(macierz, DNA, x, k)       #Przyjmowanie jako indeksow wyniku dzialania funkcji
        if (Score(Motifs, DNA, k) > Score(bestMotifs, DNA, k)):   #Jesli Score od Motifs wiekszy niz Score od BestMotifs
            bestMotifs = Motifs                                   #Motifs przyjmuje miano BestMotifs
        else:                                                  #Jesli Score od aktualnych bestMotifs jest wiekszy
            return bestMotifs                                  #Zakoncz dzialanie nieskonczonej petli i zwroc bestMotifs

if __name__ == "__main__":
    dna = ['acgtatcttc', 'tttatcggg', 'gactataaaa']
    print("\nDziałanie randomized_motif_search dla dna: ", dna)
    best = randomized_motif_search(dna, 3, len(dna))
    print("Best Motif: ", best)
    print("Score: ", Score(best, dna, 3))

    dna = ['tagtggtcttttgagtgtagatctgaagggaaagtatttccaccagttcggggtcacccagcagggcagggtgacttaat',
           'cgcgactcggcgctcacagttatcgcacgtttagaccaaaacggagttggatccgaaactggagtttaatcggagtcctt',
           'gttacttgtgagcctggttagacccgaaatataattgttggctgcatagcggagctgacatacgagtaggggaaatgcgt',
           'aacatcaggctttgattaaacaatttaagcacgtaaatccgaattgacctgatgacaatacggaacatgccggctccggg',
           'accaccggataggctgcttattaggtccaaaaggtagtatcgtaataatggctcagccatgtcaatgtgcggcattccac',
           'tagattcgaatcgatcgtgtttctccctctgtgggttaacgaggggtccgaccttgctcgcatgtgccgaacttgtaccc',
           'gaaatggttcggtgcgatatcaggccgttctcttaacttggcggtgcagatccgaacgtctctggaggggtcgtgcgcta',
           'atgtatactagacattctaacgctcgcttattggcggagaccatttgctccactacaagaggctactgtgtagatccgta',
           'ttcttacacccttctttagatccaaacctgttggcgccatcttcttttcgagtccttgtacctccatttgctctgatgac',
           'ctacctatgtaaaacaacatctactaacgtagtccggtctttcctgatctgccctaacctacaggtcgatccgaaattcg']
    print("Działanie randomized_motif_search dla dna: ", dna)
    best = randomized_motif_search(dna, 30, len(dna))
    print("Best Motif: ", best)
    print("Score: ", Score(best, dna, 30))