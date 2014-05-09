USE `bio_mcyamaha`;

CREATE TABLE `sequence`
(
	`id`
		INT
		PRIMARY KEY
		NOT NULL
		AUTO_INCREMENT,
	`year`
		INT
		NOT NULL,
	`consensus`
		BOOL
		NOT NULL
		DEFAULT false,
	`description`
		VARCHAR(32)
		NOT NULL
		DEFAULT '',
	`sequence`
		TEXT
		NOT NULL
);