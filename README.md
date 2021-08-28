## Teoria i Metody Optymalizacji
Przedmiot o nazwie 'Teoria i Metody Optymalizacji' był realizowany na studiach magisterkich na kierunku Automatyka i Robotyka, na specjalności 'Technologie Informacyjne w Systemach Automatyki' pod przewodnictwem doktor Ewy Szlachcic. Celem przedmiotu było poznanie teorii i dostępnych obecnie metod optymalizacji lokalnej oraz globalnej. Zajęcia składały się z wykładu, na którym przedstawiana była teoria oraz z projektu, którego celem było dogłębne poznanie oraz zaimplementowanie jednej z wybranych metod optymalizacji.

## Platforma
Wspierany jest jedynie system operacyjny Linux. Sama aplikacja była testowana tylko i wyłącznie na Ubuntu 20.04 LTS.

## Skład projektu
````
├── app
├── documents
├── latex
````

## Wymagania
Wymagane biblioteki:
* sympy
* numpy
* matplotlib
* math
* tkinter

## Opis projektu
Celem projektu było minimalizowanie funkcji f(x) nieliniowej, ciągłej na zbiorze ograniczeń X. Za metodę rozwiązywania została przyjęta metoda **Powella-Gaussa-Seidela**. Natomiast ze metodę minimum w kierunku została wybrana metoda **Złotego Podziału**. Metody te to algorytmy optymalizacji lokalnej. W celu lepszego zrozumienia tych metod należy rozbić je na trzy osobne tematy:

1. Metoda przesuwanej funkcji kary Powella.
2. Minimalizacja funkcji metodą Gaussa-Seidela.
3. Minimalizacja w kierunku metodą Złotego Podziału.

Aby w pełni zrozumieć działanie algorytmów zachęcam do przeczytania dokumentów z katalogu **documents**. Katalog ten zawiera:
* **Sprawozdanie z projektu** - które zawiera opis, każdej z wyżej wymienionych metod, przetestowane przykłady wraz z numerycznymi wynikami. Dodatkowo została opisana aplikacja graficzna umożliwiająca rysowanie funkcji celu wraz z ograniczeniami.
* **Skrypt** - z którego korzystano w celu zrozumienia zasady działania metod od strony matematycznej.

#### Wygląd aplikacji
![image](https://i.imgur.com/awlP3Qn.png)

#### Przykłady zastosowań
| ![Screen 1](https://i.imgur.com/GXyrszT.png) | ![Screen 2](https://i.imgur.com/GV92EIJ.png) |
|----------------------------------------------|----------------------------------------------|
| ![Screen 3](https://i.imgur.com/l1odZva.png) | ![Screen 4](https://i.imgur.com/TA3lRQ5.png) |

## Kontakt
* Autorzy: Patrycja i Kamil Kiełbasa.
* Email: kamilkielbasa64@gmail.com