USE `bio_mcyamaha`;

CREATE TABLE `sequence`
(
	`id`
		INT
		PRIMARY KEY
		NOT NULL
		AUTO_INCREMENT,
	`subtype`
		VARCHAR(16)
		NOT NULL
		DEFAULT '',
	`year`
		INT
		NOT NULL,
	`consensus`
		BOOL
		NOT NULL
		DEFAULT false,
	`threshold`
		DECIMAL(3, 2)
		DEFAULT NULL,
	`georegion`
		VARCHAR(32)
		DEFAULT NULL,
	`country`
		VARCHAR(32)
		DEFAULT NULL,
	`description`
		VARCHAR(32)
		NOT NULL
		DEFAULT '',
	`sequence`
		TEXT
		NOT NULL
);
