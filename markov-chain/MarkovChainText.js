if (! this.MarkovChainText)
{
	this.MarkovChainText = {};
}
	
(function ()
{
	this.MarkovChainText =
	{
		// public:
		
		init:
			function ()
			{
				// Train our Markov model on some source text.
				// I've chosen the first five chapters of 
				// Jane Austen's _Pride and Prejudice_.
				
				MarkovChainText.train( $('.TRAIN').val() );
				
				// Generate some fake text based on the source text.
				
				MarkovChainText.emit(MarkovChainText._LENGTH);
			},
			
		train:
			function (text)
			{
				// Partition the strings into words-like sequences.
				
				var words = text.split(/\s*(?=\b)/);
				
				for (var i = 0; i < words.length - 1; i++)
				{
					// Punctuation, numbers, and other incidentals are
					// going to generate noise that looks un-source-like,
					// but we're more interested in simplicity right now
					// than spookily meaningful text generation.
					
					var current = words[i].replace(/\"/, '');
					var next    = words[i + 1].replace(/\"/, ''); 
					
					// We're storing the strings as keys in an object
					// literal, which may not like keys consisting of
					// punctuation and quotes.  Escape it for now.
					
					current     = encodeURIComponent(current);
					next        = encodeURIComponent(next);
					
					MarkovChain.train(current, next);
				}
			},
			
		emit:
			function (length)
			{
				var current = MarkovChain.start();
				var emitted = 0;
				var text = '';
				
				while (emitted < length)
				{
					// Generate and unescape the next word generated
					// by the model based solely on the current word.
					
					text += decodeURIComponent(current) + ' ';
					emitted++;
					
					current = MarkovChain.next(current);
				}
				
				// Write the generated text to the browser.
				
				var span = 
					$('<span></span>')
						.hide()
						.attr('class', 'TEXT')
						.text(text);
						
				$('body').append(span);
				span.fadeIn('slow');
			},
			
		// private:
		
		_LENGTH: 500
	}
})();