#test packages
scDAPP::r_package_test()

#can use this to easily run pipeline to current directory
projdir = 'path/to/TestCovidDataset/'


#run pipeline with options
scDAPP::scRNAseq_pipeline_runner(
	       datadir = paste0(projdir, '/datadir/'),
	       outdir = paste0(projdir, '/outs/TESTUN/'),
               sample_metadata = paste0(projdir, '/sample_metadata.csv'),
               comps = paste0(projdir, '/comps.csv'),
               Pseudobulk_mode = T, #set to F if no replicates

               use_labeltransfer = T,
	       refdatapath = paste0(projdir, '/labeltransferref/LabelTransferRef_SCTnormalized.rds'),
	       m_reference = paste0(projdir, '/labeltransferref/LabelTransferRefMarkers.rds'),

               species = 'Homo sapiens',

               workernum = 1,
               input_seurat_obj = T #set to F for cellranger input
               )
