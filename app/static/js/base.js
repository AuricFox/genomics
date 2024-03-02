// ===================================================================
// Toggles display for sequence action descriptions
// ===================================================================
// Default to Sequence Analysis


$(document).ready(function(){
    var active = $('.active').attr('id')
    var desc = '#' + active + '-desc';

    $('.description').not($(desc)).hide();  // Hide descriptions
});

$('.seq-action').hover(function(){
    var id = $(this).attr('id');            // Get id of clicked element
    var desc = '#' + id + '-desc';

    $('.description').not($(desc)).hide();  // Hide other descriptions
    $(desc).show();                         // Show the selected description
});


