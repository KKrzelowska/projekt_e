def przyklad_score(s, DNA, k):
    """Oblicza wynik Score dla motywu
        (ilość identycznych nukleotydów na danych pozycjach)
        s - lista indeksow
        DNA - tabelka sekwencji nukleotydowych
        k - dlugosc meru"""

    score = 0                                    #ustawienie poczatkowej wartosci score
    for i in range(k):                           #iterowanie po dlugosci meru
        cnt = dict(zip("acgt", (0, 0, 0, 0)))    #przypisanie liczb nukleotydom
        for j, s_val in enumerate(s):            #iterowanie po dlugosci i elementach s
            base = DNA[j][s_val + 1]               #wyciagniecie nukleotydu z DNA
            cnt[base] += 1                       #zwiekszenie wartosci dla danego nukleotydu we wczesniej stworzonym slowniku
        score += max(cnt.values())               #dodawanie do score najwiekszej wartosci ze slownika
    return score

def Score(s, DNA, k):
    """Oblicza wynik Score dla motywu
        (ilość identycznych nukleotydów na danych pozycjach)
         s - lista indeksow
         DNA - tabelka sekwencji nukleotydowych
         k - dlugosc meru"""
    score = 0
    actg = [0, 0, 0, 0]
    for x in range(0, k):                  #iterowanie po długosci meru
        for i in range(0, len(DNA)):       #iterowanie po ilosci sekwencji
            if(DNA[i][s[i] + x] == 'a'):   #Sprawdzanie tożsamości nukleotydów
                actg[0] += 1
            if (DNA[i][s[i] + x] == 'c'):
                actg[1] += 1
            if (DNA[i][s[i] + x] == 't'):
                actg[2] += 1
            if (DNA[i][s[i] + x] == 'g'):
                actg[3] += 1
        score += max(actg)                #sumowanie największej liczby powtorzen nukleotydow na danej pozycji
        actg = [0,0,0,0]                  #czyszczenie zmiennej
    return score


def profil(motif, DNA, k):
    """"Zwraca profil prawdopodobienstwa
        motif - motyw o dlugosci k
        DNA - tabelka sekwencji nukleotydowych
        k - dlugosc meru"""
    macierz = []
    for x in range(0, k):                     #Pętla od długosci meru
        tab = [0, 0, 0, 0]                    #stworzenie kolumny dla wystepowania ACGT
        for y in range(0, len(motif)):        #Iterowanie po ilosci sekwencji w motif
            if DNA[y][x + motif[y]] == 'a':   #Jeśli w sekwencji nr y na miejscu x + indeks nukleotydem jest adenina
                tab[0] += 1                   #Wtedy wartosc kolumny adeniny zwieksza sie o 1
            if DNA[y][x + motif[y]] == 'c':
                tab[1] += 1
            if DNA[y][x + motif[y]] == 'g':
                tab[2] += 1
            if DNA[y][x + motif[y]] == 't':
                tab[3] += 1
        for i in range(0, len(tab)):            #Pętla po długości kolumny
            tab[i] = (tab[i] + 1) / len(motif)  #zastosowanie reguły Cromwella
        macierz.append(tab)
    return macierz


def podobienstwo(macierz, DNA, nr, k):
    """"Zwraca indeks najbardziej podobnego meru z sekwencji o zadanym indeksie
        macierz - wynik dzialania funkcji profil
        DNA - tabelka sekwencji nukleotydowych
        nr - numer sekwencji
        k - dlugosc meru"""
    Najpd = 0                                #Najwyzsze podobienstwo wynosi 0
    indeks = 0
    for x in range(0, len(DNA[nr]) - k):     #Pętla po dlugosci sekwencji zmniejszonej o dlugosc meru
        pd = 1                               #Podobienstwo wynosi 1
        for i in range(0, k):                #Petla po dlugosci meru
            if DNA[nr][x + i] == 'a':        #Jeśli w danej sekwencji na miejscu x + 1 nukleotydem jest adenina
                pd *= macierz[i][0]          #pomnóż podobienstwo przez wartosc prawdopodobienstwa
            if DNA[nr][x + i] == 'c':
                pd *= macierz[i][1]
            if DNA[nr][x + i] == 'g':
                pd *= macierz[i][2]
            if DNA[nr][x + i] == 't':
                pd *= macierz[i][3]
            if pd > Najpd:                   #Jesli podobienstwo jest wieksze od najwyzszego podobienstwa
                Najpd = pd                   #Nadaj Najpd nową wartosc
                indeks = x                   #zapisz indeks meru w specjalnej zmiennej
    return indeks


if __name__ == "__main__":
    dna = ['acgtatcaaa', 'aaaatcggg', 'gactataaaa']
    print("\nDla dna: ", dna)
    print("Score [0, 0, 0]: ",Score([0, 0, 0], dna, 3))
    print("Score [3, 3, 3]: ", Score([3, 3, 3], dna, 3))
    macierz = profil([1,2,3], dna, 4)
    print("profil: ", macierz)
    pd = podobienstwo (macierz, dna, 2, 4)
    print("podobienstwo: ", pd)

    dna = ['gggatgtatcttc',
           'gggtggatcgggg',
           'ggggactataaaa']

    print("\nDla dna: ", dna)
    print("Score [0, 0, 0]: ", Score([0, 0, 0], dna, 3))
    print("Score [3, 3, 3]: ", Score([3, 3, 3], dna, 3))
    macierz = profil([1, 2, 3], dna, 4)
    print("profil: ", macierz)
    pd = podobienstwo(macierz, dna, 2, 4)
    print("podobienstwo: ", pd)
