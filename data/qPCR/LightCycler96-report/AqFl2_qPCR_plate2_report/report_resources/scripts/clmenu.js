// JavaScript for collapsing / expanding menu
function toggleMenu(objID) {
	if (!document.getElementById) return;
	var ob = document.getElementById(objID).style;
	var img = document.getElementById("expandStateIcon");
	if (ob.display == 'block')
	{
		ob.display = 'none';
		img.src = "report_resources/images/expand_icon.png";
	}
	else
	{
		ob.display = 'block';
		img.src = "report_resources/images/collapse_icon.png";
	}
}

function Initialize(linkid){
var elem = document.getElementById(linkid);
if (typeof elem.onclick == "function") {
    elem.onclick.apply(elem);
}

if (elem.tagName == "A")
{
  report = document.getElementById("report_viewer");
  report.src = elem.href;
}
}

function selectItem(itemId, parentContainerId)
{
  var elem = document.getElementById(parentContainerId);
  if (elem)
  {	
     var items = elem.getElementsByTagName("a");
	 for (var i=0; i < items.length; i++) 
	 {
        if (items[i].id == itemId)
		{
			items[i].className  = "selected";
		}
		else
		{
			items[i].className  = "";
		}
	 }
  }
}