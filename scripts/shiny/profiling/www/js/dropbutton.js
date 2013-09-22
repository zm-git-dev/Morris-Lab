
// Copyright (C) 2013 by Chris Warth

// jQuery code for connecting a bootstrap dropdown button to a Shiny control.

// Event handling code to relay events to the shiny binding

// $(document).ready(function(){
//     $('.dropdown-menu').on('click', 'li', function(evt){
// 	// evt.target is the button that was clicked
// 	var el = $(evt.currentTarget);
// 	var dropButton = $(evt.delegateTarget);

// 	console.log("on click .dropdown-menu")
// 	console.log("    data() = " + $(el).data())
// 	console.log($(el))
//         $(dropButton).data('choice', $(el).data('choice'))

// 	// Raise an event to signal that the value changed
// 	dropButton.trigger("change");
//     });
// });


// Establish a Shiny binding for dropmenu button controls.

var dropButtonBinding = new Shiny.InputBinding();
console.log("extending inputbinding for dropButtonBinding")
$.extend(dropButtonBinding, {
    find: function(scope) {
	console.log("find for dropbuttonBinding")
	var tags = $(scope).find('.dropdown-menu');
	console.log(tags);
	return tags;
    },
    getId: function(el) {
	console.log("getId for dropbuttonBinding");
	console.log("    element: " + el);
	var id = (Shiny.InputBinding.prototype.getId.call(this, el) || el.name);
	console.log("    Id: " + id);
	return id;
    },


    getValue: function(el) {	
	return $(el).data("btn-value");
    },
    setValue: function(el, value) {
	console.log("setValue for dropbuttonBinding")
	$(el).val(value);
    },
    subscribe: function(el, callback) {
	console.log("subscribe for dropbuttonBinding")
	// $(el).on('change.dropButtonBinding', function(event) {
	//     console.log("on change dropButtonBinding")
	//     console.log("    event " + event)
        //     var i = $(el).data('i') || 0;
        //     $(el).data('i', i + 1);
        //     $(el).data('val', $(el).data('choice') + "-" + i);

	//     $('.dropdown.open').removeClass('open');
        //     callback();
	//     console.log("exiting change.dropButtonBinding")
	// });
	$(el).on("click.dropButtonBinding", "a[data-value]", function(e) {
	    var i = $(el).data('i') || 0;
            $(el).data('i', i + 1);
	    $(el).data("btn-value", $(e.target).data("value") + "-" + i);
	    $('.dropdown.open').removeClass('open');
	    callback();
	});
    },
    unsubscribe: function(el) {
	console.log("unsubscribe for dropbuttonBinding")
	$(el).off('.dropButtonBinding');
    }
});
Shiny.inputBindings.register(dropButtonBinding, 'shiny.dropButton');

