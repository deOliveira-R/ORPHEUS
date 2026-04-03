/**
 * Sortable tables for ORPHEUS Sphinx documentation.
 *
 * Any <table> with class "sortable" becomes click-to-sort by column header.
 * Features:
 *   - Click column header to sort (toggles ascending/descending)
 *   - Active sort column header is highlighted
 *   - Category group separators (thick borders) between distinct values
 *     in the sorted column
 */
$(document).ready(function () {
  $("table.sortable").each(function () {
    var $table = $(this);
    var $headers = $table.find("thead th");
    var numCols = $headers.length;

    $headers.css("cursor", "pointer").attr("title", "Click to sort");

    function applySortVisuals(colIdx) {
      /* Insert thick border between rows where the sorted column value changes. */
      var $tbody = $table.find("tbody");
      var $rows = $tbody.find("tr");

      // Clear previous group borders
      $rows.removeClass("group-first");

      var prevVal = null;
      $rows.each(function () {
        var val = $(this).children("td").eq(colIdx).text().trim();
        if (prevVal !== null && val !== prevVal) {
          $(this).addClass("group-first");
        }
        prevVal = val;
      });
    }

    $headers.click(function () {
      var $th = $(this);
      var colIdx = $th.index();
      var $tbody = $table.find("tbody");
      var rows = $tbody.find("tr").toArray();
      var ascending = !$th.data("sort-asc");

      // Sort rows by the text content of the target column
      rows.sort(function (a, b) {
        var aText = $(a).children("td").eq(colIdx).text().trim();
        var bText = $(b).children("td").eq(colIdx).text().trim();

        // Try numeric sort first
        var aNum = parseFloat(aText);
        var bNum = parseFloat(bText);
        if (!isNaN(aNum) && !isNaN(bNum)) {
          return ascending ? aNum - bNum : bNum - aNum;
        }
        // Fall back to string sort
        return ascending
          ? aText.localeCompare(bText)
          : bText.localeCompare(aText);
      });

      $tbody.empty().append(rows);

      // Update sort indicator and active-column highlight
      $headers
        .removeData("sort-asc")
        .removeClass("sort-active")
        .find(".sort-indicator")
        .remove();
      $th.data("sort-asc", ascending).addClass("sort-active");
      $th.append(
        '<span class="sort-indicator"> ' +
          (ascending ? "\u25B2" : "\u25BC") +
          "</span>"
      );

      // Apply group separator borders
      applySortVisuals(colIdx);
    });

    // Trigger default sort on first column (Method) on page load
    $headers.eq(1).trigger("click");
  });
});
