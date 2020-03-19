from przydatne_funkcje import Score, profil, podobienstwo
from itertools import product    #Załączenie funkcji pozwalającej iterować wszystkie mozliwości

def brute_force_motif_search(dna, t, n, l):
    """Funkcja zwracająca najczęstszy motyw po przeszukaniu wszystkich możliwości
       DNA - tabelka sekwencji nukleotydowych
       t - ilosc sekwencji w dna
       n - ilosc nukleotydow w sekwencji
       l - dlugosc meru"""
    # Danej implementacji algorytmu użyłam, gdyż jest zgodna z założeniami:
    # a)iteruje po wszystkich możliwościach
    # b)przyrównuje wynik Score do bestScore
    # c)zwraca bestMotif jako liste indeksów o najwyższym Score
    #biblioteka itertools byla wspomniana na laboratoriach
    #podaję to rozwiązanie jako alternatywę dla mój_brute.py

    bestMotif = []
    for j in range(0, t):
        bestMotif.append(0)
    bestScore = 0
    lista = []
    for i in range(0, n - l + 1):      #stworzenie listy elementów na których funkcja będzie iterować
        lista.append(i)
    kombinacje = list(product(lista, repeat = t))   #funkcja tworzy liste wszystkich kombinacji podanych elementow
    for komb in kombinacje:                         #Iterujemy po wszystkich mozliwościach
        if(Score(komb, dna, l) > bestScore):        #Porównujemy wyniki Score do bestScore
            bestScore = Score(komb, dna, l)
            bestMotif = komb
    return bestMotif

if __name__ == "__main__":
    dna = ['accggg', 'ataggg', 'aacggg', 'attggg']
    print("\nDziałanie brute_force_motif_search dla dna: ", dna)
    print("Dla dlugosci meru: ", 3)
    best = brute_force_motif_search(dna, len(dna), len(dna[1]), 3)
    print("Best Motif: ", best)
    print("Score: ", Score(best, dna, 3))
    print("Najlepszy motyw: ", dna[0][best[0]:best[0] + 3])


    dna = ['gggatgtatcttc',
           'gggtggatcgggg',
           'ggggactataaaa']
    print("\nDziałanie brute_force_motif_search dla dna: ", dna)
    print("Dla dlugosci meru: ", 3)
    best = brute_force_motif_search(dna, len(dna), len(dna[1]), 3)
    print("Best Motif: ", best)
    print("Score: ", Score(best, dna, 3))
    print("Najlepszy motyw: ", dna[0][best[0]:best[0] + 3])

