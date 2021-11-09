function createSeqInput() {
    var newDiv = document.createElement("div");
    newDiv.className = 'seqgroup';

    var tgc = document.createElement("seq");
    tgc.innerHTML="Forbidden Sequences(in IUPAC code):<input type='text' name='seq'> ";
    newDiv.appendChild(tgc);

    var btnElem = document.createElement("button");
    btnElem.type = "button";
    btnElem.textContent = "delete";
    btnElem.addEventListener("click", removeUrlBox);
    newDiv.appendChild(btnElem);

    var element = document.getElementById("moreseq");
    element.appendChild(newDiv);
}

function removeUrlBox() {
    this.closest('.seqgroup').remove();
}
