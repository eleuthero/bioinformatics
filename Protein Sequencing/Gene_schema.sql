USE 'bio_mcyamaha';

CREATE TABLE 'protein'
(
	'id'
		INT
		PRIMARY KEY
		NOT NULL
		AUTO_INCREMENT,
	'fragment'
		VARCHAR(32)
		NOT NULL,
	'start'
		INT
		DEFAULT NULL,
	'stop'
		INT
		DEFAULT NULL
);