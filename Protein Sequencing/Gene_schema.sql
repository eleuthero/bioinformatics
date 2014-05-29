USE 'bio_mcyamaha';

CREATE TABLE 'protein'
(
	'id'
		INT
		PRIMARY KEY
		NOT NULL
		AUTO_INCREMENT,
	'feature'
		VARCHAR(32)
		NOT NULL,
	'start'
		INT
		DEFAULT NULL,
	'stop'
		INT
		DEFAULT NULL,
    'type'
        VARCHAR(32)
        DEFAULT NULL,
    'category'
        VARCHAR(32)
        DEFAULT NULL,
    'parent region'
        VARCHAR(32)
        DEFAULT NULL,
    'reference'
        VARCHAR(32)
        DEFAULT NULL,
    'note'
        VARCHAR(32)
        DEFAULT NULL
);