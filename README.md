# Monte Carlo metodų taikymas statistinių hipotezių tikrinimui genų raiškos duomenyse

**Autorius:** Margiris Antanas Malakauskas  
**Vadovas:** lekt. Irus Grinis  
**Vilniaus universitetas, Matematikos ir informatikos fakultetas, 2026**

---

## Projekto aprašas

Šis projektas lygina du statistinius metodus diferencinės genų raiškos analizėje:
- **Welch t-testas** — klasikinis parametrinis metodas
- **Monte Carlo permutacinis testas** — neparametrinis metodas pagrįstas atsitiktiniu modeliavimu

Analizei naudojami viešai prieinami Alzheimerio ligos hipokampo genų raiškos duomenys
(GSE1297) iš NCBI GEO duomenų bazės. Lyginami 9 sveikų kontrolinių asmenų ir 7 sunkios
Alzheimerio stadijos pacientų mėginiai, iš viso 22 283 genai.

---

## Reikalavimai

- **R** (4.0 ar naujesnė versija)
- **RStudio** (rekomenduojama)
- R paketai: 


---

## Failų struktūra

```
├── Kursinis_darbas.R
├── README.md
```

---

## Kodo struktūra

Kursinis_darbas.R buvo suskirstytas į 12 sekcijų, kurios kiekvienos reikšmė ir atliekamos funkcijos surašytos komentaruose:

Sekcijos:
1. Paketų įkėlimas
2. Duomenų atsiuntimas ir patikrinimas
3. Kontrolinių ir sunkių AD samplų išskyrimas
4. Neapdorotų duomenų peržiūra
5. Log2 transformacija
6. Vieno geno analizė (Genas 1: 1007_s_at)
7. Visų 22 283 genų analizė
8. Daugybinio tikrinimo korekcija (Benjamini-Hochberg FDR)
9. Metodų palyginimas
10. Specifinio (tarp meodų labiausiai nesutarto) geno analizė: Genas 203894_at
11. Skaičiavimo laiko palyginimas
12. Permutacijų skaičiaus stabilumo analizė

---

## Paleidimas

Atidaryti `Kursinis_darbas.R` RStudio aplinkoje ir paleisti visą kodą iš viršaus į apačią.

> Dėl Monte Carlo metodo, kodas gali ilgai užrukti (daugiau nei valandą)

---

## Pagrindiniai rezultatai

- Spearman koreliacija tarp metodų: **ρ = 0.9948**
- Rekomenduojamas permutacijų skaičius: **B = 5 000**
- Monte Carlo greitis: ~**1 150 kartų lėtesnis** už t-testą

---

## Šaltiniai

- Blalock et al. (2004) — GSE1297 duomenų rinkinys
- Benjamini & Hochberg (1995) — FDR korekcija
- Phipson & Smyth (2010) — permutacinių p-reikšmių korekcija
- Metropolis & Ulam (1949) — Monte Carlo metodas
