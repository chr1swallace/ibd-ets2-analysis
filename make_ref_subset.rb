#!/home/cew54/localc/bin/ruby
require 'zlib'
require 'multiple_files_gzip_reader'
require 'tempfile'
require_relative 'dirs.rb'

block=ARGV[0] # chr*_block*
infile=ARGV[1] # list of chr pos
outfile=ARGV[2] # plink root

## get and check list of snps
lines=File.readlines(infile)
chr=nil
lines.each { |l|
  fields=l.chomp.split("\t")
  if !chr.nil? && chr!=fields[0]
    abort "multiple chromosomes seen in #{infile}"
  end
  chr=fields[0]
}

vcffile="#{DIR}/reference/byblock/#{block}.vcf.gz"
puts "subsetting vcf file #{vcffile}, searching for #{lines.length} SNPs"

samplefile="/home/cew54/share/Data/reference/1000GP_Phase3/sparse_basis/EUR.sample" # for now, can make an option later

vcftemp= `mktemp /tmp/tempvcf.XXXXXX`.chomp
command = "zcat #{vcffile} | " +
          "sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | " +
          "#{ENV['HOME']}/localc/bin/vcftools " +
          " --gzvcf - " +
          " --remove-indels --recode --remove-filtered-all --keep #{samplefile} " +
          " --positions #{infile} --stdout > #{vcftemp} "
system(command)

## next: make plink
command = "#{ENV['HOME']}/localc/bin/plink --vcf #{vcftemp} --out #{outfile}"
system(command)

## cleanup
File.unlink(vcftemp)

exit 0;
