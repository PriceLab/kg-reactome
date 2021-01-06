vizmap = [

   {selector:"node", css: {
       "text-valign":"top",
       "text-halign":"center",
       "background-color": "lightgreen",
       "border-color": "black",
       "content": "data(label)",
       "border-width": "1px"
   }},

    {selector:"node[type='molecule']", css: {
       "text-valign":"center",
       "text-halign":"center",
	"shape": "ellipse",
        "width": "120px;",
	"height": "120px;",
	"background-color": "#D0FFD0"
        }},

    {selector:"node[type='metabolite']", css: {
	"shape": "ellipse",
        "width": "120px;",
	"height": "80px;",
	"background-color": "lightblue"
        }},

    {selector:"node[type='drug']", css: {
	"shape": "round-heptagon",
        "width": "80px;",
	"height": "80px;",
	"background-color": "#D0D0D0",
	"border-width": "2px",
	"border-color": "black"
        }},

    {selector:"node[type='reaction']", css: {
       "text-valign":"center",
       "text-halign":"center",
	"shape": "round-diamond",
        "width": "200px;",
	"height": "100px;",
	"background-color": "lightblue"
        }},

    {selector:"node[type='protein']", css: {
	"shape": "ellipse",
        "width": "50px;",
	"height": "50px;",
	"background-color": "orange",
	"font-size": "20px"
        }},

    {selector:"node[type='complex']", css: {
	"shape": "barrel",
        "width": "100px;",
	"height": "100px;",
	"background-color": "beige",
	"font-size": "20px"
        }},

    {selector:"node[type='pathway']", css: {
	"shape": "round-triangle",
        "width": "60px;",
	"height": "60px;",
	"background-color": "lightblue",
	"font-size": "10px"
        }},

    {selector:"node[type='complexMember']", css: {
        "text-valign":"center",
        "text-halign":"center",
	"shape": "ellipse",
        "width": "30px;",
	"height": "30px;",
	"background-color": "white",
	"font-size": "8px"
        }},

   {selector:"node:selected", css: {
       "text-valign":"center",
       "text-halign":"center",
       "border-color": "black",
       "font-size": "32px",
       "content": "data(id)",
       "border-width": "3px",
       "overlay-opacity": 0.5,
       "overlay-color": "blue"
       }},


   {selector: 'edge', css: {
       'line-color': 'maroon',
       'target-arrow-shape': 'triangle',
       'target-arrow-color': 'black',
       'curve-style': 'bezier',
       'label': "data(interaction)"
       }},
    
   {selector: "edge[interaction='modifies']", css: {
       'line-color': 'gray',
       'arrow-scale': 2,
       'target-arrow-shape': 'triangle',
       'target-arrow-color': 'black',
       'line-style': "dashed",
       'curve-style': 'bezier'
   }},

   {selector: "edge[interaction='substrateOf'][reversible='true']", css: {
       'line-color': 'red',
       'source-arrow-shape': 'triangle',
       'source-arrow-color': 'black',
       'arrow-scale': 2,
       'target-arrow-shape': 'triangle',
       'target-arrow-color': 'black',
       'line-style': "solid",
       'curve-style': 'bezier'
       }},

   {selector: "edge[interaction='substrateOf'][reversible='false']", css: {
       'line-color': 'red',
       'arrow-scale': 2,
       'target-arrow-shape': 'triangle',
       'target-arrow-color': 'black',
       'line-style': "solid",
       'curve-style': 'bezier'
       }},


   {selector: "edge[interaction='produces'][reversible='true']", css: {
       'line-color': 'green',
       'arrow-scale': 2,
       'target-arrow-shape': 'triangle',
       'target-arrow-color': 'black',
       'source-arrow-shape': 'triangle',
       'source-arrow-color': 'black',
       'line-style': "solid",
       'curve-style': 'bezier'
       }},

   {selector: "edge[interaction='produces'][reversible='false']", css: {
       'line-color': 'green',
       'arrow-scale': 2,
       'target-arrow-shape': 'triangle',
       'target-arrow-color': 'black',
       'line-style': "solid",
       'curve-style': 'bezier'
       }}
    

]

