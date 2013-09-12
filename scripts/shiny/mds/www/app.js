
// Copyright (C) 2013 by Chris Warth

// jQuery code for connecting a bootstrap dropdown button to a Shiny control.


// Event handling code to relay events to the shiny binding

$(document).ready(function(){
    $('.dropdown-menu').on('click', 'li', function(evt){
	// evt.target is the button that was clicked
	var el = $(evt.target);

	console.log(".dropdown-menu click")
	console.log(evt.target.id)
	// Raise an event to signal that the value changed
	el.trigger("change");
    });
});


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
	console.log(el);
	var id = (Shiny.InputBinding.prototype.getId.call(this, el) || el.name);
	console.log(id);
	return id;
    },
    getValue: function(el) {
	console.log("getValue for dropbuttonBinding")
	return $(el).val();
    },
    setValue: function(el, value) {
	console.log("setValue for dropbuttonBinding")
	$(el).val(value);
    },
    getState: function(el) {
	// Store options in an array of objects, each with with value and label
	var options = new Array(el.length);
	console.log("getState for dropbuttonBinding")
	for (var i = 0; i < el.length; i++) {
            options[i] = { value:    el[i].value,
			   label:    el[i].label,
			   selected: el[i].selected };
	}

	return {
            label: $(el).parent().find('label[for=' + el.id + ']').text(),
            value:    this.getValue(el),               
            options:  options
	};
    },
    receiveMessage: function(el, data) {
	var $el = $(el);

	console.log("receiveMessage for dropbuttonBinding")
    },
    subscribe: function(el, callback) {
	console.log("subscribe for dropbuttonBinding")
	$(el).on('change.dropButtonBinding', function(event) {
	    console.log("callback dropButtonBinding")
            callback();
	});
    },
    unsubscribe: function(el) {
	console.log("unsubscribe for dropbuttonBinding")
	$(el).off('.dropButtonBinding');
    }
});
Shiny.inputBindings.register(dropButtonBinding, 'shiny.dropButton');

