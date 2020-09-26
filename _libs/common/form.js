/**
 * Return an object of parameter values
 * 
 * @param elements      Elements as extracted using document.querySelectorAll('form')[0]
 * @return              ObjectOfInputs
 */
var formDataExtractor = elements => [].reduce.call(elements, (objectOfInputs, element) => {
    // If it is a number, but in string form, convert it to a number
    if (!isNaN(element.value)) {
        objectOfInputs[element.name] = parseFloat(element.value);
    } else {
        objectOfInputs[element.name] = element.value;
    }
    return objectOfInputs;
}, {});

/**
 * Read inputs from the form and enter them into an object.
 * 
 * @params              None.
 * @return              An object containing all the form inputs.
 */
function readInputs() {
    // Form elements
    var form = document.querySelectorAll('form')[0];
    var objectOfInputs = formDataExtractor(form);
    return objectOfInputs;
}