/*  burkhac@students.wwu.edu
 *  CSCI 474 -- Bioinformatics
 *  No rights reserved.  Share and enjoy !
 *
 *  This javascript library is used to demonstrate
 *  set covers over SNPs.
 */

if (! this.Snp)
{
	this.Snp = {};
}

(function ()
{
	this.Snp =
	{
		// public:

		DEFAULT_PATTERNCOUNT: 7,		
		
		DEFAULT_SNPCOUNT: 9,
		
		init:
			function ()
			{
				$('input:text.SCOUNT')
					.val(Snp.DEFAULT_SNPCOUNT)
					.on
					(
						'change',
						function ()
						{
							Snp._buildTable();
						}
					);
					
				$('input:text.PCOUNT')
					.val(Snp.DEFAULT_PATTERNCOUNT)
					.on
					(
						'change',
						function ()
						{
							Snp._buildTable();
						}
					);

				$('input:button.RESET')
					.on
					(
						'click',
						function ()
						{
							Snp._buildTable();
						}
					);
					
				$('input:button.COVER')
					.on
					(
						'click',
						function ()
						{
							Snp._setCoverGreedy();
						}
					);
				
				Snp._buildTable();
								
				$('table.SNP')
					.on
					(
						'change',
						'input:checkbox',
						Snp._onTagChange
					);
			},
	
		// private:

		_distinguishability: [],
		
		_snpCount:
			function /* int */ ()
			{
				return parseInt($('input:text.SCOUNT').val()) || Snp.DEFAULT_SNPCOUNT;
			},

		_patternCount:
			function /* int */ ()
			{			
				return parseInt($('input:text.PCOUNT').val(), 10) || Snp.DEFAULT_PATTERNCOUNT;
			},
			
		_buildTable:
			function ()
			{
				$('table.SNP')
					.empty();
				
				Snp._buildRowHeader();
				
				for (var i = 1; i <= Snp._snpCount(); i++)
				{
					Snp._buildRow(i);
				}
				
				Snp._distinguishability = [];
			},

		_buildRowHeader:
			function ()
			{
				var tr =
					$('<tr></tr>')
						.append( $('<td></td>') );
						
				for (var i = 1; i <= Snp._patternCount(); i++)
				{
					$('<th></th>')
						.html('P<sub>' + i + '</sub>')
						.appendTo(tr);
				}
				
				// Column for checkboxes.
				
				tr.append( $('<th></th>') );

				$('table.SNP').append(tr);
			},
			
		_buildRow:
			function (i)
			{
				var tr =
					$('<tr></tr>')
						.append
						(
							$('<td class="LEFT"></td>')
								.html('S<sub>' + i + '</sub>')
						);
						
				for (var i = 0; i < Snp._patternCount(); i++)
				{
					var sClass =
						Math.random() > 0.5
							? "MAJOR"
							: "MINOR";
							
					$('<td></td>')
						.append
						(
							$('<div></div>')
								.addClass(sClass)
						)
						.appendTo(tr);
				}
				
				tr.append
				(
					$('<td></td>')
						.append( $('<input type="checkbox"></input') )
				)
				
				$('table.SNP').append(tr);
			},
			
		_onTagChange:
			function (e)
			{		
				$('table.SNP tr')
					.removeClass('SELECTED');

				$('input:checkbox:checked')
					.closest('tr')
						.addClass('SELECTED');
						
				$('.DISTINGUISHABLE').removeClass('DISTINGUISHABLE');
				$('.INDISTINGUISHABLE').removeClass('INDISTINGUISHABLE');

				// Determine whether the checked snp collection
				// is a cover.
				
				Snp._isCover();
			},
			
		// Returns true if the checked collection of SNPs is a cover.
		
		_isCover:
			function ()
			{
				// Generate integer representation of elements of
				// P for the selected candidate subset of S.
				
				var aP = [];
				
				for (var j = 0; j < Snp._patternCount(); j++)
					aP[j] = 0;
				
				// Over the candidate cover ...
				
				var tags =
					$('input:checkbox:checked')
						.closest('tr');
				
				for (var i = 0; i < tags.length; i++)
				{
					for (var j = 0; j < Snp._patternCount(); j++)
					{
						var div = tags.eq(i).find('div:eq(' + j + ')');
						
						if (div.hasClass('MAJOR'))
						{
							aP[j] = (aP[j] | 0) + (1 << i);
						}
					}
				}
				
				// console.log('integer representations: ' + aP);
								
				// We now have an integer representation of the profiles.
				// Can we distinguish them all ?  
				// Look for clashes in the array; those correspond
				// to non-distinguishable patterns.  Also, note that
				// the algorithm I'm using is stupid and is O(n^2).
				// We could do better, but this is just a toy.
				// Replace clashes with false values; replace unique
				// values with true values.
				
				var isCover = true;
				
				for (var i = 0; i < aP.length; i++)
				{
					var t = aP[i];

					if ('boolean' == typeof(t))
					{
						// Already processed this element.
						
						continue;
					}
					
					// Assume ith element is distinguishable until
					// proven otherwise.
					
					aP[i] = true;

					for (var j = i + 1; j < aP.length; j++)
					{
						if (aP[j] == t)
						{
							// Proven otherwise.
							
							aP[j] = aP[i] = isCover = false;
						}
					}
				}
				
				// At this point, aP[i] is true if the ith individual
				// is distinguishable, and false if not.
				
				console.log('distinguishability: ' + aP);
				
				// Save the distinguishability matrix.  We may need
				// it again later.
				
				Snp._distinguishability = aP;
				
				// Update the UI.
				
				$.each
				(
					aP,
					function (i, b)
					{
						tags
							.find('div:eq(' + i + ')')
								.toggleClass('DISTINGUISHABLE',    b)
								.toggleClass('INDISTINGUISHABLE', !b);
								
						$('th')
							.eq(i)
								.toggleClass('DISTINGUISHABLE', b);
					}
				);
				
				return isCover;
			},
			
		_setCoverGreedy:
			function ()
			{
				if (Snp._isCover())
				{
					alert('Selected SNPs already form a cover.');
					return;
				}
				
				var aS = [];
				
				// Select our candidate S[i] from the collection of unchecked
				// SNP rows.
				
				var trs =
					$('input:checkbox:not(:checked)')
						.closest('tr');
						
				// For each S[i], determine whether it can distinguish
				// patterns P[j] from P[k] for all (j != k) in the pattern
				// set.  For each pair (j, k) that S[i] can distinguish,
				// it gets a point.  The algorithm is O(m * n^2) just for
				// one pass.  Nothing to be proud of, but again, this is just
				// a toy.
				
				for (var i = 0; i < Snp._snpCount(); i++)
				{
					aS[i] = 0;					
					var tr = trs.eq(i);
					
					for (var j = 0; j < Snp._patternCount(); j++)
					{
						for (var k = j + 1; k < Snp._patternCount(); k++)
						{
							if (Snp._distinguishability[j] || Snp._distinguishability[k])
							{
								continue;
							}
					
							var b1 = tr.find('div:eq(' + j + ')').hasClass('MAJOR');
							var b2 = tr.find('div:eq(' + k + ')').hasClass('MAJOR');
							
							// If the two candidates aren't the same flavor
							// of allele (either both major or both minor)
							// then S[i] can distinguish them.  One point
							// for you.
							
							if (b1 ^ b2)
							{
								aS[i] = (aS[i] || 0) + 1;
							}
						}
					}
				}
				
				console.log('pairwise patterns distinguished: ' + aS);
								
				// At this point, aS[i] is the number of patterns that S[i]
				// can distinguish.  In accordance with the greedy algorithm,
				// we want to select the candidate snp with the greatest ability
				// to distinguish patterns.  Identify the first maximal value
				// in the array and select it as our candidate.
				
				var imax = 0;
				var vmax = 0;
				
				for (var i = 0; i < aS.length; i++)
				{
					if (aS[i] > vmax)
					{
						vmax = aS[i];
						imax = i;
					}
				}
				
				// We now have the index of an optimally discriminating row.
				// Select it.
				
				$('input:checkbox:not(:checked)')
					.eq(imax)
						.click();							
			}
	}
})();
