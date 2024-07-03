# interaction samples
```bash
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 19 15:40:43 ~
mv SchunterC_RNASeq_CPOS-221130-CS-15577b CS_RNASeq
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 19 15:36:09 ~/CS_RNAseq/FastQC
for i in *.zip;do unzip ${i};done
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 19 15:39:32 ~/CS_RNAseq/primary_seq
ll *.fastq.gz|perl -alne '@a=split /\_/,$F[-1];$c=$a[0]."_R".$a[1];print "$F[-1]\t$c"' >../data_list.txt
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 19 15:42:40 ~/CS_RNAseq
find -name "fastqc_data.txt"|perl -ne 'chomp;print "$_\t"'|perl -alne 'print "cat $_ > all.samples.fastqc_data.txt"' >1.sh
sh 1.sh
less Overrep_seq.txt|perl -alne '$i++;$nm=Over.$i;print ">$nm\n$_"' > Overrep_seq.txt.1
mv Overrep_seq.txt.1 Overrep_seq.txt

less data_list.txt|cut -f 2|perl -alne '@a=split /\./, $_;$b=$a[0].".fq.gz";print "$_\t$b" ' >data_list.txt.1
mv data_list.txt.1 data_list.txt

# (base) kang1234@celia-PowerEdge-T640 Thu Oct 19 15:50:12 ~/CS_RNAseq
cp ~/New_caledonia/quality_control.pl ./
mv FastQC fastqc1

vi quality_control.pl # this script still needs to be revised to rename the raw fastq file
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 19 16:06:23 ~/CS_RNAseq
nohup perl quality_control.pl --input data_list.txt --raw_dir primary_seq --trim_dir ~/software/Trimmomatic-0.39 --kraken_lib ~/software/kraken2/library --fastqc ~/software/FastQC/fastqc >quality_control.process 2>&1 &
# 5458

# de novo
# jlkang@hpc2021 Fri Oct 27 15:20:09 /lustre1/g/sbs_schunter/Kang
mkdir CS_RNAseq
nohup scp Abu*_1.fq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/CS_RNAseq > nohup_Abu1.out 2>&1
nohup scp Abu*_2.fq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/CS_RNAseq > nohup_Abu2.out 2>&1
nohup scp M*_1.fq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/CS_RNAseq > nohup_M1.out 2>&1
nohup scp M*_2.fq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/CS_RNAseq > nohup_M2.out 2>&1
```

```Trinity.cmd
#!/bin/bash
#SBATCH --job-name=Run_trinity      # 1.    Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlkang@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=28
#SBATCH --mem=250G                    # 6.    Request total amount of RAM
#SBATCH --time=7-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=%x_%j.err             #    Standard error log as $job_name_$job_id.err

module load trinity/2.13.2
module load bowtie2/2.4.4
module load rsem/1.3.1_test
Trinity --seqType fq --full_cleanup --min_kmer_cov 1 --left Abu11_1.fq.gz,Abu12_1.fq.gz,Abu14_1.fq.gz,Abu16_1.fq.gz,Abu17_1.fq.gz,Abu18_1.fq.gz,Abu20_1.fq.gz,Abu21_1.fq.gz,Abu22_1.fq.gz,Abu24_1.fq.gz,Abu26_1.fq.gz,Abu28_1.fq.gz,Abu29_1.fq.gz,Abu31_1.fq.gz,Abu3_1.fq.gz,Abu32_1.fq.gz,Abu35_1.fq.gz,Abu38_1.fq.gz,Abu39_1.fq.gz,Abu40_1.fq.gz,Abu42_1.fq.gz,Abu43_1.fq.gz,Abu44_1.fq.gz,Abu45_1.fq.gz,Abu46_1.fq.gz,Abu49_1.fq.gz,Abu50_1.fq.gz,Abu51_1.fq.gz,Abu5_1.fq.gz,Abu52_1.fq.gz,Abu54_1.fq.gz,Abu56_1.fq.gz,Abu59_1.fq.gz,Abu60_1.fq.gz,Abu6_1.fq.gz,Abu8_1.fq.gz --right Abu11_2.fq.gz,Abu12_2.fq.gz,Abu14_2.fq.gz,Abu16_2.fq.gz,Abu17_2.fq.gz,Abu18_2.fq.gz,Abu20_2.fq.gz,Abu21_2.fq.gz,Abu22_2.fq.gz,Abu24_2.fq.gz,Abu26_2.fq.gz,Abu28_2.fq.gz,Abu29_2.fq.gz,Abu31_2.fq.gz,Abu32_2.fq.gz,Abu3_2.fq.gz,Abu35_2.fq.gz,Abu38_2.fq.gz,Abu39_2.fq.gz,Abu40_2.fq.gz,Abu42_2.fq.gz,Abu43_2.fq.gz,Abu44_2.fq.gz,Abu45_2.fq.gz,Abu46_2.fq.gz,Abu49_2.fq.gz,Abu50_2.fq.gz,Abu51_2.fq.gz,Abu52_2.fq.gz,Abu5_2.fq.gz,Abu54_2.fq.gz,Abu56_2.fq.gz,Abu59_2.fq.gz,Abu60_2.fq.gz,Abu6_2.fq.gz,Abu8_2.fq.gz --CPU 28 --max_memory 240G --bflyCPU 4 --bflyHeapSpaceMax 20G --output trinity_Abu 
Trinity --seqType fq --full_cleanup --min_kmer_cov 1 --left M108_1.fq.gz,M111_1.fq.gz,M114_1.fq.gz,M122_1.fq.gz,M123_1.fq.gz,M127_1.fq.gz,M13_1.fq.gz,M133_1.fq.gz,M134_1.fq.gz,M151_1.fq.gz,M156_1.fq.gz,M165_1.fq.gz,M166_1.fq.gz,M173_1.fq.gz,M176_1.fq.gz,M18_1.fq.gz,M23_1.fq.gz,M25_1.fq.gz,M26_1.fq.gz,M27_1.fq.gz,M28_1.fq.gz,M29_1.fq.gz,M30_1.fq.gz,M34_1.fq.gz,M35_1.fq.gz,M36_1.fq.gz,M40_1.fq.gz,M45_1.fq.gz,M56_1.fq.gz,M63_1.fq.gz,M64_1.fq.gz,M68_1.fq.gz,M7_1.fq.gz,M75_1.fq.gz,M76_1.fq.gz,M78_1.fq.gz,M83_1.fq.gz,M84_1.fq.gz,M86_1.fq.gz,M88_1.fq.gz,M91_1.fq.gz,M92_1.fq.gz,M93_1.fq.gz,M94_1.fq.gz,M95_1.fq.gz,M96_1.fq.gz,M97_1.fq.gz --right M108_2.fq.gz,M111_2.fq.gz,M114_2.fq.gz,M122_2.fq.gz,M123_2.fq.gz,M127_2.fq.gz,M13_2.fq.gz,M133_2.fq.gz,M134_2.fq.gz,M151_2.fq.gz,M156_2.fq.gz,M165_2.fq.gz,M166_2.fq.gz,M173_2.fq.gz,M176_2.fq.gz,M18_2.fq.gz,M23_2.fq.gz,M25_2.fq.gz,M26_2.fq.gz,M27_2.fq.gz,M28_2.fq.gz,M29_2.fq.gz,M30_2.fq.gz,M34_2.fq.gz,M35_2.fq.gz,M36_2.fq.gz,M40_2.fq.gz,M45_2.fq.gz,M56_2.fq.gz,M63_2.fq.gz,M64_2.fq.gz,M68_2.fq.gz,M7_2.fq.gz,M75_2.fq.gz,M76_2.fq.gz,M78_2.fq.gz,M83_2.fq.gz,M84_2.fq.gz,M86_2.fq.gz,M88_2.fq.gz,M91_2.fq.gz,M92_2.fq.gz,M93_2.fq.gz,M94_2.fq.gz,M95_2.fq.gz,M96_2.fq.gz,M97_2.fq.gz --CPU 28 --max_memory 240G --bflyCPU 4 --bflyHeapSpaceMax 20G --output trinity_M
```

```bash
# jlkang@hpc2021 Fri Oct 27 16:21:04 /lustre1/g/sbs_schunter/Kang/CS_RNAseq
sbatch Trinity.cmd
```

```bash
# de novo in SNORLAX
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 30 09:29:15 ~/CS_RNAseq/kraken
cat Abu*_1.fq.gz > AbuT_1.fq.gz
cat Abu*_2.fq.gz > AbuT_2.fq.gz
export SHARED_DIR=$PWD
function docker_drap_cmd { echo "docker run --rm -e LOCAL_USER_ID=`id -u $USER` -u 1003:1003 -v $SHARED_DIR:$SHARED_DIR -w `pwd` sigenae/drap "; }
`docker_drap_cmd`runDrap -o Abu_runDrap -R1 AbuT_1.fq.gz -R2 AbuT_2.fq.gz -s FR --dbg trinity --dbg-mem 400 --no-trim --norm-mem 400 --no-rate --write
nohup `docker_drap_cmd`runDrap -o Abu_runDrap \
          -1 AbuT_1.fq.gz \
          -2 AbuT_2.fq.gz \
          -s FR --dbg trinity \
          --dbg-mem 400 --no-trim --norm-mem 400 --no-rate --run > Abu_runDrap.process 2>&1 &
# [1] 20444

cat M*_1.fq.gz > MTotal_1.fq.gz
cat M*_2.fq.gz > MTotal_2.fq.gz
export SHARED_DIR=$PWD
function docker_drap_cmd { echo "docker run --rm -e LOCAL_USER_ID=`id -u $USER` -u 1003:1003 -v $SHARED_DIR:$SHARED_DIR -w `pwd` sigenae/drap "; }
`docker_drap_cmd`runDrap -o M_runDrap -R1 MTotal_1.fq.gz -R2 MTotal_2.fq.gz -s FR --dbg trinity --dbg-mem 400 --no-trim --norm-mem 400 --no-rate --write
nohup `docker_drap_cmd`runDrap -o M_runDrap \
          -1 MTotal_1.fq.gz \
          -2 MTotal_2.fq.gz \
          -s FR --dbg trinity \
          --dbg-mem 400 --no-trim --norm-mem 400 --no-rate --run > M_runDrap.process 2>&1 &
# [1] 9334
```

```script_trinity.cmd
#!/bin/bash
#SBATCH --job-name=Trinity      # 1.    Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2.    Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlkang@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd           # 3.    Request a partition
#SBATCH --qos=long                 # 4.    Request a QoS
#SBATCH --ntasks=1                  # 5.    Request total number of tasks (MPI workers)
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G                     # 6.    Request total amount of RAM
#SBATCH --time=14-00:00:00             # 7.    Job execution duration limit day-hour:min:sec
#SBATCH --output=%x_%j.out            # 8.    Standard output log as $job_name_$job_id.out
#SBATCH --error=%x_%j.err             #    Standard error log as $job_name_$job_id.err


cat *_1.fq.gz > all_1.fq.gz
cat *_2.fq.gz > all_2.fq.gz
module load trinity/2.13.2
module load bowtie2/2.4.4
module load rsem/1.3.1_test
module load samtools/1.14
Trinity --seqType fq --full_cleanup --min_kmer_cov 1 --left all_1.fq.gz --right all_2.fq.gz --CPU 32 --max_memory 240G --bflyCPU 4 --bflyHeapSpaceMax 20G
```

```bash
# jlkang@hpc2021 Wed Apr 10 15:06:28 /lustre1/g/sbs_schunter/Kang/JD_data/Thermal_rates/kraken
sbatch script_trinity.cmd
```

```bash
# orthologous detection
# (base) kang1234@celia-PowerEdge-T640 Fri Jun 28 15:05:15 ~/CS_RNAseq/kraken
mkdir orthologue
cp M_runDrap/e-rmbt_editing/all_contigs.second_pass.fa orthologue/M_assembly.fa
cp Abu_runDrap/e-rmbt_editing/all_contigs.second_pass.fa orthologue/Abu_assembly.fa
# 然后用TransDecoder得到其蛋白质序列
nohup TransDecoder.LongOrfs -t Abu_assembly.fa  > Abu_transdecoder.process 2>&1 &
# [1] 25157
nohup TransDecoder.LongOrfs -t M_assembly.fa > M_transdecoder.process 2>&1 &
# [2] 25320

# 把整理好的pep序列放入一个文件夹
# (base) kang1234@celia-PowerEdge-T640 Fri Jun 28 15:59:59 ~/CS_RNAseq/kraken/orthologue/orthofinder_input_pep
mkdir orthofinder_input_pep/
cd orthofinder_input_pep/
less ../Abu_assembly.fa.transdecoder_dir/longest_orfs.pep|perl -alne 'if (/>/){$i++;my $nm=Abu."_"."$i";print">$nm"}else{print"$_"}' > Abu_pep.fasta
less ../M_assembly.fa.transdecoder_dir/longest_orfs.pep|perl -alne 'if (/>/){$i++;my $nm=M."_"."$i";print">$nm"}else{print"$_"}' > M_pep.fasta
# (base) kang1234@celia-PowerEdge-T640 Fri Jun 28 16:03:34 ~/CS_RNAseq/kraken/orthologue
nohup orthofinder -f orthofinder_input_pep -a 30 >orthofinder-process 2>&1 &
# [1] 5796
```
