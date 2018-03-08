# Алгоритм Нидлмана — Вунша 
Был предложен в 1970 году Солом Нидлманом и Кристианом Вуншем для выравнивания двух аминокислотных или нуклеотидных последовательностей.
[Подробнее](https://ru.wikipedia.org/wiki/%D0%90%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC_%D0%9D%D0%B8%D0%B4%D0%BB%D0%BC%D0%B0%D0%BD%D0%B0_%E2%80%94_%D0%92%D1%83%D0%BD%D1%88%D0%B0)
##  Использование
* '-m', '--mistake', Penalty for gap or insert/Штраф за ошибку, введите значение float, по умолчанию значение = -5
* '-W', '--W_output', Print W, Вывести матрицу
* '-l', '--label', The Needleman–Wunsch algorithm to align protein sequences and build heat-map.
* '-i', '--input', Input filename for analise/Имя выходного файла, обязательно

## Пример 
$ python3 Dz3-arparse.py -i input_file.vcf -m -7 -W -l

Последовательности из файла input_file.vcf будут выравнены со штрафом -7, выведена матрица выравниваний и её карта-плотности(?)

$ python3 Dz3-arparse.py -i input_file.vcf -W

Последовательности из файла input_file.vcf будут выравнены со стандартным штрафом -5, выведена только матрица выравниваний 
