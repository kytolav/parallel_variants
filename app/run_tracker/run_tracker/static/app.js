$(document).ready(function() {

    function ajaxQuery(queryData, successCallback, path) {
      // Sending AJAX request to server
      $.ajax({
        url: path,
        type: "POST",
        dataType   : 'json',
        contentType: 'application/json; charset=UTF-8',
        data: JSON.stringify(queryData),
        success: function(result) {
            console.log('Successful ajax');
            successCallback(result);
        }
      });
    };

    // Initializing the table
    $('#run-table').DataTable( {
        "data": {}
    } );

    // Update table function
    function updateTable(data) {
        var datatable = $('#run-table').dataTable().api();
        datatable.clear();
        datatable.rows.add(data);
        datatable.draw();
        console.log(data);
    };

    // Periodical AJAX for fetching the data table contents
    var intervalID = setInterval(function(){
        // Fetching new table data from server
        ajaxQuery({'get_data': 1}, updateTable, "http://compbio.uta.fi:8008/table") //"http://localhost:8008/table")
    }, 5000);
} );