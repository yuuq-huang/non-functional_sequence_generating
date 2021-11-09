function createInput() {
    var newDiv = document.createElement("div");
    newDiv.className = 'group';

    var tgc = document.createElement("tgc");
    tgc.innerHTML="Temp:<input type='' name='temp'> GC-content:<input type=' ' name='GC'> range:<input type='' name='ini'> - <input type=' ' name='end'>";
    newDiv.appendChild(tgc);

    var btnElem = document.createElement("button");
    btnElem.type = "button";
    btnElem.textContent = "delete";
    btnElem.addEventListener("click", removeUrlBox);
    newDiv.appendChild(btnElem);

    var element = document.getElementById("tempgc");
    element.appendChild(newDiv);
}

function removeUrlBox() {
    this.closest('.group').remove();
}
function add() {
    var tgc = document.createElement('tgc')
    tgc.innerHTML="<div id='divline'> <p>Temp:<input type='' name='temp'> GC-content:<input type=' ' name='GC'> range:<input type='' name='ini'> - <input type=' ' name='end'> <input type='button' value='delete' onclick='deline('divline'+a)'> </p> </div>";
    document.getElementById("tempgc").appendChild(tgc)
}
function deline(divid) {
    var li = document.getElementById(divid);
    li.parentNode.removeChild(li);

}