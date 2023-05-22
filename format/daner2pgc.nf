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
		.map {it -> [(it.name =~ /basic\.(.+)\.(MDD.+)_(.+)\.num\.xls/)[0][1..3], it]}
		
	DANER_CH = Channel
		.fromPath(params.daner)
		.map {it -> [(it.name =~ /daner_(.+)\.(MDD.+)_(.+)\.gz/)[0][1..3], it]}
		
	TEMPLATE_CH = Channel
		.fromPath(params.template)
		
	CFF_CH = Channel
		.fromPath(params.cff)
		
	FAI_CH = Channel
		.fromPath(params.fai)
	
	DANER_BASIC_CH =
	BASIC_CH
		.join(DANER_CH)
		
	DANER_PGC_CH = DANER_BASIC_CH
		.combine(TEMPLATE_CH)
		.combine(CFF_CH)
		.combine(FAI_CH)
		
	PGC_CH = PGC(DANER_PGC_CH)
}

process PGC {
	tag "${analysis.join(' ')}"
	
	publishDir "results", mode: 'copy'
	
	cpus = 1
	memory = 4.GB
	time = '10m'
	
	input:
	tuple val(analysis), path(basic), path(daner), path(template), path(cff), path(fai)
	
	output:
	path("*.tsv.gz")
	
	script:
	"""
	pgc.R --name="MDD symptom: ${analysis[1]} (${analysis[2]}) [${analysis[0]}]" \
	--daner=${daner} \
	--basic=${basic} \
	--analyst="${params.analyst}" \
	--fai=${fai} \
	--cff=${cff} \
	--template=${template} \
	--out=pgc-${analysis[1]}-${analysis[2]}-${analysis[0]} \
	--pop=EUR
	
	gzip *.tsv
	"""
	
}