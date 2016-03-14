var homzScale = d3.scale.linear()
    .domain( [0,         .33,        .66,        1])
    .range(['#66bd63', '#fee08b', '#f46d43', '#d73027'])


var freqScale = d3.scale.linear()
.domain( [0,         .33,        .66,        1])
.range(['#66bd63', '#fee08b', '#f46d43', '#d73027'])


var domainGen = function(ar){
var start = ar[0]
var end = ar[1]
var range = end - start
return [start, start + (0.33*range), start + (0.66*range), end]
}

//Waits to execute the javascript code 200 miliseconds. This is needed because shiny loads the javascript before it loads the table
//and thus has nothing to color.

freqTableVals = [] // initialize array to store table values for color scaling.
homzTableVals = []

window.setInterval(function() {
freqTableVals = []
homzTableVals = []
//Grab the allele freq column's max and mins
  d3.selectAll('#asf_table tbody tr td:nth-child(4)')
.each(function() {
var cellValue = d3.select(this).text();
freqTableVals = _.union(freqTableVals, [parseFloat(cellValue)])
})

freqColorRange = d3.extent(freqTableVals)
freqScale.domain(domainGen(freqColorRange))

//grab the homz column's max and mins
  d3.selectAll('#asf_table tbody tr td:nth-child(5)')
.each(function() {
var cellValue = d3.select(this).text();
homzTableVals = _.union(homzTableVals, [parseFloat(cellValue)])
})

homzColorRange = d3.extent(homzTableVals)
homzScale.domain(domainGen(homzColorRange))

  //Set the colors for the allele freq
  d3.selectAll('#asf_table tbody tr td:nth-child(4)')
    .style('background-color', function() {
      var cellValue = d3.select(this).text();
      return (freqScale(cellValue))
    })
  //set colors for homz freq
  d3.selectAll('#asf_table tbody tr td:nth-child(5)')
    .style('background-color', function() {
      var cellValue = d3.select(this).text();
      return (homzScale(cellValue))
    })
}, 200);
