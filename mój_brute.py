from przydatne_funkcje import Score, profil, podobienstwo

def m_brute_force_motif_search(DNA, t, n, l):
    """Funkcja zwracająca najczęstszy motyw po przeszukaniu wszystkich możliwości
       DNA - tabelka sekwencji nukleotydowych
       t - ilosc sekwencji w dna
       n - ilosc nukleotydow w sekwencji
       l - dlugosc meru"""
    #Danej implementacji algorytmu użyłam, gdyż jest zgodna z założeniami:
    # a)iteruje po wszystkich możliwościach
    # b)przyrównuje wynik Score do bestScore
    # c)zwraca bestMotif jako liste indeksów o najwyższym Score
    #Zasada działania jest taka sama

    bestMotif = []
    for j in range(0, t):
        bestMotif.append(0)
    bestScore = 0             #początkowo bestScore o jak najmniejszej wartosci
    s = []
    for x in range(0, t):     #pętla iterująca przez ilosc sekwencji w DNA
        s.append(0)           #zapełnianie s początkowymi indeksami ustalonymi na 0
    for y in range(0, (pow((n-l)+1, t))):               #Pętla iterująca przez wszystkie mozliwości (usunełam wstawione poprzednio błędnie -1)
        if(Score(s, DNA, l) > bestScore):       #Porownywanie wyniku Score do bestScore
            for i in range(0, len(s)):
                bestMotif[i] = s[i]
            bestScore = Score(s, DNA, l)
        if(int(s[0]) % (n - l) == 0 and s[0] != 0): #Jesli reszta z dzielenia pierwszego elementu s(zmieniającego się najczęściej) przez długoścSekwencji/dlugoscMotywu równa 0, zakaz dzielenia prez zero
            for z in range(1, t):                   #Pętla od 1 do ilosci sekwencji w dna
                if s[z] == n - l:                   #Jeśli indeks danej sekwencji = długoścSekwencji - dlugoscMotywu
                    s[z] = 0                        #indeks zerujemy
                else:                               #Jeśli indeks jeszcze nie osiągnął maksymalnej wartości
                    s[z] += 1                       #Zwiększ go o 1
                    break                           #i wyjdz z petli
            s[0] = 0                                #Wyzeruj pierwszy indeks
        else:                                       #Jesli indeks pierwszego elementu jeszcze nie osiągnął maksymalnej wartości
            s[0] += 1                               #Zwiększ go o jeden
    return bestMotif

if __name__ == "__main__":
    dna = ['accggg', 'atgggg', 'aacggg', 'attggg']
    print("\nDziałanie brute_force_motif_search dla dna: ", dna)
    print("Dla dlugosci meru: ", 3)
    best = m_brute_force_motif_search(dna, len(dna), len(dna[1]), 3)
    print("Best Motif: ", best)
    print("Score: ", Score(best, dna, 3))
    print("Najlepszy motyw: ", dna[0][best[0]:best[0] + 3])


    dna = ['gggatgtatcttc',
           'gggtggatcgggg',
           'ggggactataaaa']
    print("\nDziałanie brute_force_motif_search dla dna: ", dna)
    print("Dla dlugosci meru: ", 3)
    best = m_brute_force_motif_search(dna, len(dna), len(dna[1]), 3)
    print("Best Motif: ", best)
    print("Score: ", Score(best, dna, 3))
    print("Najlepszy motyw: ", dna[0][best[0]:best[0]+3])

