USE 'bio_mcyamaha';

CREATE TABLE 'Proteins'
        (
            'HXB2_K03455 '
				VARCHAR(3)
				NOT NULL
				DEFAULT '',
            'HXB2_Position'
				VARCHAR(10)
				NOT NULL
				DEFAULT '',
			'frame1' 
				VARCHAR(5)
				DEFAULT NULL,
            'frame2' 
				VARCHAR(5)
				DEFAULT NULL,
            'frame3' 
				VARCHAR(5)
				DEFAULT NULL,
            'aminoAcid' 
				VARCHAR(5)
				DEFAULT NULL,
            'aminoAcidDes' 
				VARCHAR(50)
				DEFAULT NULL,
            'proteinNumber' 
				VARCHAR(20)
				DEFAULT NULL,
            'gene' 
				VARCHAR(30)
				DEFAULT NULL,
            'protein' 
				VARCHAR(70)
				DEFAULT NULL,
            'RNAFeature' 
				VARCHAR(100)
				DEFAULT NULL,
            'proteinFeature' 
				VARCHAR(70)
				DEFAULT NULL,
            'secondFrame' 
				VARCHAR(5) 
				DEFAULT NULL,
            'numbering' 
				VARCHAR(10)
				DEFAULT NULL,
            'proteinFeature2' 
				VARCHAR(50)
				DEFAULT NULL,
            'thirdFrame' 
				VARCHAR(10)
				DEFAULT NULL,
            'numbering3' 
				VARCHAR(6)
				DEFAULT NULL,
            'proteinFeature3' 
				VARCHAR(70)
				DEFAULT NULL,
            'reference' 
				VARCHAR(50)
				DEFAULT NULL,
		);
