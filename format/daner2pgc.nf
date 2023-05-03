/* Convert daner sumstats to PGC VCF-like rich header format */
nextflow.enable.dsl=2

params.basic = "../meta/distribution/**/basic.*.num.xls"
params.daner = "../meta/distribution/**/daner_*.gz"
params.cff = "../CITATION.cff"
params.analyst = "Mark James Adams, UoE"
params.template = "pgc.glue"
params.fai = "http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.fai"

workflow {
	BASIC_CH = Channel
		.fromPath(params.basic)
		.map {it -> [(it.name =~ /basic\.(.+)\.num\.xls/)[0][1], it]}
		
	DANER_CH = Channel
		.fromPath(params.daner)
		.map {it -> [(it.name =~ /daner_(.+)\.gz/)[0][1], it]}
		
	TEMPLATE_CH = Channel
		.fromPath(params.template)
		
	CFF_CH = Channel
		.fromPath(params.cff)
		
	FAI_CH = Chanenl
		.fromPath(params.fai)
	
	DANER_BASIC_CH =
	BASIC_CH
		.join(DANER_CH)
		
	PGC_CH = PGC(DANER_BASIC_CH, TEMPLATE_CH, CFF_CH, FAI_CH)
}

process {
	tag "{analysis}"
	
	cpus = 1
	memory = 4.GB
	time = '10m'
	
	input:
	each tuple val(analysis), path(basic), path(daner)
	path(template)
	path(cff)
	path(fai)
	
	output:
	path("*.tsv.gz")
	
	script:
	"""
	Rscript pgc.R --name=MDD1 \
	--daner=${daner} \
	--basic=${basic} \
	--analyst="${params.analyst}" \
	--fai={fai} \
	--cff={cff} \
	--template={glue} \
	--out=pgc-mdd1-comm \
	--name="MDD symptom: low mood (MDD1)" \
	--pop=EUR
	"""
	
}