// ===================================================================
// Toggles display for sequence action descriptions
// ===================================================================
// Default to Sequence Analysis
$('.description').not($('#analysis-desc')).hide();  // Hide descriptions
$('.forms').not($('#analysis-form')).hide();        // Hide forms

$('.seq-action').click(function () {
    var id = $(this).attr('id');            // Get id of clicked element
    
    var desc = '#' + id + '-desc';
    var form = '#' + id + '-form';
    
    $('.description').not($(desc)).hide();  // Hide other descriptions
    $(desc).show();                         // Show the selected description

    $('.forms').not($(form)).hide();        // Hide other forms
    $(form).show();                         // Show the selected form
});