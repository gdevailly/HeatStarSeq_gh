$(function() {
	$(document).on({
		'shiny:disconnected': function(event) {
			alert('Connection lost, sorry :-(\n\n' +
				  'Click OK to come back to Heat*seq main page\n\n' +
				  'If you see this message again, please try relaunching the app after freeing some resources on your computer, or try with a more powerful computer.\n' +
				  'Contact: guillaume.devailly _at_ rolsin.ed.ac.uk'
			);
			window.location = 'http://www.heatstarseq.roslin.ed.ac.uk/';
		}
	});
});


