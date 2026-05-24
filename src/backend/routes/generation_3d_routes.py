from flask import Blueprint, jsonify, request
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from dotenv import load_dotenv
import base64, asyncio, io, traceback
from PIL import Image
from playwright.async_api import async_playwright

load_dotenv()

generation_3d_routes = Blueprint("generation_3d_routes", __name__)

RENDER_W     = 800
RENDER_H     = 500
RENDER_SCALE = 2

def normalize_to_inchi(input_str: str) -> str:
    if not input_str:
        raise ValueError("Empty input")
    
    input_str = input_str.strip()
    
    if input_str.startswith("InChI="):
        return input_str
    
    try:
        mol = Chem.MolFromSmiles(input_str)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {input_str}")
        
        inchi = Chem.MolToInchi(mol)
        if not inchi:
            raise ValueError(f"Could not convert SMILES to InChI: {input_str}")
        
        return inchi
    
    except Exception as e:
        raise ValueError(f"Invalid input (not SMILES or InChI): {input_str} - {e}")

def inchi_to_mol3d(inchi):
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        return None, None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        return None, None
    AllChem.UFFOptimizeMolecule(mol)
    return mol, Chem.MolToMolBlock(mol)


def build_html(sdf, vw, vh):
    return (
        "<!DOCTYPE html><html><head><style>"
        "* { margin:0; padding:0; box-sizing:border-box; }"
        "body { background:white; width:" + str(vw) + "px; height:" + str(vh) + "px; overflow:hidden; }"
        "#v { width:" + str(vw) + "px; height:" + str(vh) + "px; }"
        "</style></head><body>"
        "<div id='v'></div>"
        "<script src='https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js'></script>"
        "<script>"
        "var v = $3Dmol.createViewer('v', { backgroundColor: 'white', antialias: true });"
        "v.addModel(" + repr(sdf) + ", 'sdf');"
        "v.setStyle({}, { stick: { radius: 0.12, colorscheme: 'Jmol' }, sphere: { scale: 0.22, colorscheme: 'Jmol' } });"
        "v.zoomTo();"
        "v.render();"
        "window.__done = true;"
        "</script></body></html>"
    )


async def playwright_render(sdf):
    vw = RENDER_W * RENDER_SCALE
    vh = RENDER_H * RENDER_SCALE
    html = build_html(sdf, vw, vh)

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        ctx = await browser.new_context(
            viewport={"width": vw, "height": vh},
            device_scale_factor=RENDER_SCALE
        )
        page = await ctx.new_page()
        await page.set_content(html, wait_until="networkidle")
        await page.wait_for_function("window.__done === true", timeout=15000)

        try:
            await page.wait_for_function(
                """() => {
                    var c = document.querySelector('canvas');
                    if (!c || c.width === 0 || c.height === 0) return false;
                    var ctx2d = c.getContext('2d');
                    if (!ctx2d) return false;
                    var d = ctx2d.getImageData(0, 0, Math.min(c.width,100), Math.min(c.height,100)).data;
                    for (var i = 0; i < d.length; i += 4) {
                        if (d[i] < 245 || d[i+1] < 245 || d[i+2] < 245) return true;
                    }
                    return false;
                }""",
                timeout=6000
            )
        except Exception:
            await page.wait_for_timeout(800)

        png = await page.screenshot(
            type="png",
            clip={"x": 0, "y": 0, "width": vw, "height": vh}
        )
        await browser.close()

    img = Image.open(io.BytesIO(png)).resize((RENDER_W, RENDER_H), Image.LANCZOS)

    pixels = list(img.getdata())
    non_white = sum(1 for px in pixels if px[0] < 245 or px[1] < 245 or px[2] < 245)
    if non_white < 100:
        raise ValueError("Rendered image is blank — molecule did not render in time")

    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


def run_async(coro):
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    try:
        return loop.run_until_complete(coro)
    finally:
        loop.close()


@generation_3d_routes.route("/api/generate_3d", methods=["POST"])
def generate_3d():
    try:
        input_str = request.json.get("inchi")  
        
        inchi = normalize_to_inchi(input_str)
        
        mol, sdf = inchi_to_mol3d(inchi)
        if mol is None:
            return jsonify({"error": "Could not generate 3D molecule"}), 400
        
        return jsonify({"sdf": sdf, "inchi": inchi})  
    
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@generation_3d_routes.route("/api/render_3d_image", methods=["POST"])
def render_3d_image():
    try:
        input_str = request.json.get("inchi")
        
        inchi = normalize_to_inchi(input_str)
        
        mol, sdf = inchi_to_mol3d(inchi)
        if mol is None:
            return jsonify({"error": "Could not generate 3D molecule"}), 400

        try:
            png = run_async(playwright_render(sdf))
            b64 = base64.b64encode(png).decode("utf-8")
            return jsonify({"image": b64, "sdf": sdf, "inchi": inchi})

        except Exception as e:
            print("Playwright render failed:", e)
            traceback.print_exc()
            try:
                mol2d = Chem.RemoveHs(mol)
                AllChem.Compute2DCoords(mol2d)
                drawer = rdMolDraw2D.MolDraw2DCairo(RENDER_W, RENDER_H)
                drawer.DrawMolecule(mol2d)
                drawer.FinishDrawing()
                b64 = base64.b64encode(drawer.GetDrawingText()).decode("utf-8")
                return jsonify({"image": b64, "sdf": sdf, "inchi": inchi})
            
            except Exception as e2:
                return jsonify({"error": str(e2)}), 500
                
    except Exception as e:
        return jsonify({"error": str(e)}), 500