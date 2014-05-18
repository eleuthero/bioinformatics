if (! this.MarkovChain)
{
	this.MarkovChain = {};
}
	
(function ()
{
	this.MarkovChain =
	{
		train:
		
			// Iteratively train the Markov model by giving it an element
			// and an element in the source sequence immediately following it.
			// You will want to call this repeatedly, passing in the ith and
			// the (i + 1)th element for every 0 <= i < (n - 1) in the
			// source sequence.
			
			function (current, next)
			{
				// Get the transition states for the current element.
				
				var transition = MarkovChain._TRANSITIONS[current];
				
				if (undefined === transition)
				{
					// If there are none, we haven't seen this element yet.
					// Set up a transition for it.
				
					transition = {};
					transition[next] = 0;
					MarkovChain._TRANSITIONS[current] = transition;
					MarkovChain._ELEMENTS[MarkovChain._ELEMENTS.length] = current;
				}

				// Ok.  Now that we're sure we have a transition, increment the
				// number of times current transitions to next.
				
				transition[next] = (transition[next] || 0) + 1;
			},
		
		start:
		
			// Start the Markov chain.
			
			function ()
			{
				// Choose a random element to begin with ...
				
				var r = Math.floor( Math.random() * MarkovChain._ELEMENTS.length );
				return MarkovChain._ELEMENTS[r];
			},
		
		next:
		
			// This is the heart of the Markov chain.  Based on its training,
			// it generates the next element in the sequence based only on the
			// transitions seen in the source sequence.
			
			function (current)
			{
				var choices = [];
			
				// Get the transition record for the current element.
				
				var transition = MarkovChain._TRANSITIONS[current];
				
				// What elements does it transition to ?
				
				for (var choice in transition)
				{
					// How many times does current transition to this element ?
					
					var count = transition[choice];
					
					// The probability that current will now transition to this
					// element is proportional to the number of times the source
					// sequence transitions from current to this choice.
					
					// Note: this for loop is baloney.  There are better ways of
					// weighting the average than an O(N^2) algorithm, but we
					// are just keeping the code simple for now.
					
					for (var i = 0; i < count; i++)
					{
						choices[choices.length] = choice;
					}
				}

				// Ok, now we have a table of possible transition candidates.
				
				var result;
				
				if (choices.length > 0)
				{
					// Choose a choice at random since we've already accounted for their
					// frequency.
					
					result = choices[ Math.floor( Math.random() * choices.length ) ];
				}
				else
				{
					// There are some words that are dead-ends; if we run into one,
					// just restart the Markov chain.
					
					result = MarkovChain.start();
				}
	
				// if (undefined === result) debugger;
				
				return result;
			},
			
		// private:
		
		_ELEMENTS: [ ],
		_TRANSITIONS: { },
			
		_dumpTransitions:
			function ()
			{
				for (var current in MarkovChain._TRANSITIONS)
				{
					console.log(current);
					
					var transition = MarkovChain._TRANSITIONS[current];
					
					for (var next in transition)
					{
						console.log('\t' + next + ' x ' + transition[next]);
					}
				}
			}
	}
})();