if (! this.MarkovChainWord)
{
	this.MarkovChainWord = {};
}
	
(function ()
{
	this.MarkovChainWord =
	{
		// public:
		
		init:
			function ()
			{
				// Train our Markov model on some source text.
				// I've chosen the first five chapters of 
				// Jane Austen's _Pride and Prejudice_.
				
				MarkovChainWord.train( $('.TRAIN').val() );

				// Make some fake words based on the source text.
				
				MarkovChainWord.emitWords(MarkovChainWord._WORDCOUNT);
			},

		train:
			function (text)
			{
				// Remove all non-word characters and multiple whitespaces.
				
				var letters =
					text.replace(/\s+/g, ' ').replace(/[^A-Za-z ]/g, '').split('');
				
				for (var i = 0; i < letters.length - 1; i++)
				{
					var current = letters[i];
					var next    = letters[i + 1]; 
					
					MarkovChain.train(current, next);
				}
			},
		
		emitWords:
			function (count)
			{
				for (var i = 0; i < count; i++)
				{
					MarkovChainWord.emitWord();
				}
			},
			
		emitWord:
			function ()
			{
				var valid = false;
				
				// Let's ensure that our words have a minimum length
				// so that the model doesn't generate uninteresting runty
				// strings like "q".
				
				while (! valid)
				{
					var current = MarkovChain.start();
					var word = current;
				
					// Let's keep on adding to our fake word until the
					// model generates a space.
					
					while (' ' != current)
					{
						current = MarkovChain.next(current);
						word += current;
					}
					
					// Is the fake word long enough to be interesting ?
					
					valid = (word.length > MarkovChainWord._MINWORDLENGTH);
				}
				
				// Make up a random color for the fake word.
				// I couldn't do better than Soubok did it:
				// http://stackoverflow.com/questions/1152024/best-way-to-generate-a-random-color-in-javascript
				
				var color = '#' + (0x1000000 + Math.random() * 0xffffff)
					.toString(16).substr(1,6);
				
				var span = 
					$('<span></span>')
						.hide()
						.css('color', color)
						.attr('class', 'WORD')
						.text(word.toLowerCase());
						
				$('body').append(span);
				span.fadeIn('slow');
			},
			
		// private:
		
		_WORDCOUNT: 100,
		
		_MINWORDLENGTH: 4
	}
})();