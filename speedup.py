import sys

argv = sys.argv


def calc_media(lst):
    return sum(lst)/len(lst)


size = int(argv[1])
n_threads = int(argv[2])
n_repeticoes = int(argv[3])

media_seq_time = calc_media([float(x) for x in argv[4:n_repeticoes+4]])
media_par_time = calc_media([float(x) for x in argv[n_repeticoes+4:]])

print(f'Tamanho da matriz: {size}x{size}')
print(f'Média de tempo sequencial: {media_seq_time}')
print(f'Média de tempo paralelo  : {media_par_time}')
print(f'Speedup:                 : {media_seq_time/media_par_time}')
print("="*65)
