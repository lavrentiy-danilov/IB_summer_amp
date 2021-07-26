from Bio import SeqIO
import os


path = '../Data/CoralTBase/'

'''
num1 = 0
num2 = 0
os.makedirs(path + 'step1', exist_ok=True)
for f in os.listdir(path + 'cds_database/'):
  with open(path + 'step1/' + f, 'w') as f_out:
    print(path + 'step1/' + f)
    for seq_record in SeqIO.parse(path + 'cds_database/' + f, 'fasta'):
      num1 += 1
      seq_record.seq = seq_record.seq.translate()
      length = len(seq_record.seq)
      entry = seq_record.seq.count('K') + seq_record.seq.count('R')
      ratio = entry / length
      if ratio > 0.05:
        num2 += 1
        SeqIO.write(seq_record, f_out, 'fasta')
print(num1)
print(num2)


os.makedirs(path + 'step2/', exist_ok=True)
for f in os.listdir(path + 'step1/'):
  s = 'pepstats -sequence ' + path + 'step1/' + f + ' -outfile ' + path + 'step2/' + f.split('.')[0] + '.pepstats'
  print(s)
  os.system(s)


os.makedirs(path + 'step3/', exist_ok=True)
for f in os.listdir(path + 'step2/'):
  with open(path + 'step2/' + f, 'r') as f_in:
    print(path + 'step2/' + f)
    with open(path + 'step3/' + f + 'isoel', 'w') as f_out:
      for s in f_in.readlines():
        if 'Isoelectric Point' in s:
          isoel_p = -1
          try:
            isoel_p = float(s.split(' ')[-1])
          except:
            isoel_p = -1
          f_out.write(str(isoel_p) + '\n')


num1 = 0
os.makedirs(path + 'step4/', exist_ok=True)
for f in os.listdir(path + 'step3/'):
  with open(path + 'step3/' + f, 'r') as f_peps:
    print(path + 'step3/' + f)
    isoel_arr = f_peps.readlines()
    dir = path + 'step4/' + f.split('.')[0][:20] + '.t' + str(0).zfill(3)
    print(dir)
    f_out = open(dir, 'w')
    num_of_lines = 0
    for i, seq_record in enumerate(SeqIO.parse(path + 'step1/' + f.split('.')[0] + '.fasta', 'fasta')):
      isoel_p = float(isoel_arr[i])
      if isoel_p > 12:
        num_of_lines += 1
        if not num_of_lines % 1000:
          dir = path + 'step4/' + f.split('.')[0][:20] + '.t' + str(num_of_lines // 1000).zfill(3)
          print(dir)
          f_out.close()
          f_out = open(dir, 'w')
        num1 += 1
        id = str(seq_record.id)[:25]
        seq = str(seq_record.seq).replace('*', '')
        f_out.write(id + ' N N 7 298 0.1 ' + seq + '\n')
    f_out.close()
print(num1)


os.makedirs(path + 'step5/', exist_ok=True)
for f in os.listdir(path + 'step4/'):
  s = './tango_x86_64_release -inputfile=' + path + 'step4/' + f + ' > ' + path + 'step5/' + f + 'tango'
  print(s)
  os.system(s)
  os.system('rm TR* _agg*')

'''
num1 = 0
os.makedirs(path + 'step6/', exist_ok=True)
with open(path + 'step6/final_dataset.fasta', 'w') as f_out:
  for f in os.listdir(path + 'step5/'):
    with open(path + 'step5/' + f, 'r') as f_in:
      print(path + 'step5/' + f)
      id = ''
      seq = ''
      for s in f_in.readlines():
        s = s.split(' ')
        if s[0] == 'Ficher':
          id = f.split('.')[0] + '_' + s[-1]
        elif s[0] == 'aminoacido':
          seq = s[-1]
        elif s[0] == 'AGG':
          if float(s[1]) <= 500 and float(s[-5]) <= 25 and 75 <= float(s[-1]) <= 100:
            num1 += 1
            f_out.write('>' + id + seq)
            id = ''
            seq = ''
print(num1)


os.makedirs(path + 'final/', exist_ok=True)
s = 'blastp -db ../Data/database_blast -query ' + path + 'step6/final_dataset.fasta -outfmt 6 -out ' + path + 'final/final_dataset_blastp_qcov_70 -qcov_hsp_perc 70 -num_threads 5'
print(s)
os.system(s)


with open(path + 'final/final_dataset_blastp_qcov_70', 'r') as f1:
  with open(path + 'step6/final_dataset.fasta', 'r') as f2:
    fil = f2.read()
    for s in f1.readlines():
      s = s.split()[0]
      ptr1 = fil.find(s)
      ptr2 = fil.find('>', ptr1+1)
      fil = fil[:ptr1] + fil[ptr2:]
    with open(path + 'final/diff_final_dataset.fasta', 'w') as f3:
      f3.write(fil)
            
